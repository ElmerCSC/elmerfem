/*****************************************************************************/
/*
 *  Elmer, A Finite Element Software for Multiphysical Problems
 *
 *  Copyright 1st April 1995 - , CSC - IT Center for Science Ltd., Finland
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library (in file ../LGPL-2.1); if not, write
 * to the Free Software Foundation, Inc., 51 Franklin Street,
 * Fifth Floor, Boston, MA  02110-1301  USA
 *
 *****************************************************************************/

/*
 * Pin the numeric locale to "C" and, on Windows, neutralize libgfortran's
 * per-I/O locale switching.
 *
 * Background: for every formatted READ/WRITE, libgfortran temporarily forces
 * LC_NUMERIC to "C" via the C runtime's setlocale(), then restores it. On
 * targets with POSIX per-thread locales (glibc: uselocale) this is thread
 * safe. The mingw/UCRT runtime has no uselocale, so libgfortran falls back to
 * the process-global setlocale(). Under OpenMP that makes concurrent formatted
 * I/O corrupt the shared UCRT locale state, producing an intermittent crash
 * inside ucrtbase setlocale()/mbstowcs_s() (observed as a 0xC0000409 / SIGSEGV
 * during a formatted WRITE). See also the locale-related notes in
 * GeneralUtils.F90 and the MATC/Lua locale handling.
 *
 * Elmer keeps LC_NUMERIC = "C" throughout, so libgfortran's save/set-"C"/
 * restore dance is redundant here. On Windows we therefore repoint the
 * libgfortran runtime DLL's imported setlocale to a no-op that reports "C":
 * numeric formatting stays correct (it already used "C") and the unsafe global
 * locale churn disappears. Only libgfortran's import is patched; Elmer's own
 * and MATC's setlocale() calls use a separate import table and are untouched.
 *
 * The libgfortran DLL is located by enumerating the loaded modules and
 * matching any "libgfortran*.dll" (case-insensitive), rather than a hardcoded
 * soname. libgfortran has been libgfortran-5.dll since GCC 8, and GCC 16 still
 * ships version 5, but a future soname bump (libgfortran-6.dll, ...) would
 * otherwise silently disable this fix and let the intermittent crash return
 * with no diagnostic. Matching by prefix keeps the fix working across such a
 * bump. Should upstream libgfortran ever migrate mingw to a per-thread locale
 * (no setlocale import), the module is found but no thunk is patched -- a
 * harmless no-op, as the crash mechanism would no longer exist.
 *
 * Called once from the ElmerSolver library entry point, on the main thread,
 * before any OpenMP region -> the plain static guard needs no locking.
 */

#include <locale.h>

#ifdef _WIN32
#include <windows.h>
#include <string.h>
#ifdef ELMER_LOCALE_DEBUG
#include <stdio.h>
#endif

/* Stable storage for the reported locale name. */
static char elmer_c_locale[] = "C";

/* No-op replacement for libgfortran's setlocale(): never touches the real
 * runtime locale, always reports "C". */
static char *elmer_win_c_setlocale(int category, const char *locale)
{
  (void) category;
  (void) locale;
  return elmer_c_locale;
}

/* Walk one module's import table and repoint every "setlocale" thunk to
 * elmer_win_c_setlocale. Returns the number of thunks patched. */
static int elmer_patch_module_setlocale(HMODULE mod)
{
  unsigned char *base;
  IMAGE_DOS_HEADER *dos;
  IMAGE_NT_HEADERS *nt;
  IMAGE_DATA_DIRECTORY dir;
  IMAGE_IMPORT_DESCRIPTOR *desc;
  int patched = 0;

  if (!mod) return 0;
  base = (unsigned char *) mod;
  dos = (IMAGE_DOS_HEADER *) base;
  if (dos->e_magic != IMAGE_DOS_SIGNATURE) return 0;
  nt = (IMAGE_NT_HEADERS *) (base + dos->e_lfanew);
  if (nt->Signature != IMAGE_NT_SIGNATURE) return 0;
  dir = nt->OptionalHeader.DataDirectory[IMAGE_DIRECTORY_ENTRY_IMPORT];
  if (!dir.VirtualAddress) return 0;

  for (desc = (IMAGE_IMPORT_DESCRIPTOR *) (base + dir.VirtualAddress);
       desc->Name; desc++) {
    DWORD names_rva = desc->OriginalFirstThunk ? desc->OriginalFirstThunk
                                               : desc->FirstThunk;
    IMAGE_THUNK_DATA *nthunk = (IMAGE_THUNK_DATA *) (base + names_rva);
    IMAGE_THUNK_DATA *ithunk = (IMAGE_THUNK_DATA *) (base + desc->FirstThunk);

    for (; nthunk->u1.AddressOfData; nthunk++, ithunk++) {
      IMAGE_IMPORT_BY_NAME *imp;
      if (nthunk->u1.Ordinal & IMAGE_ORDINAL_FLAG) continue;
      imp = (IMAGE_IMPORT_BY_NAME *) (base + nthunk->u1.AddressOfData);
      if (strcmp((const char *) imp->Name, "setlocale") == 0) {
        DWORD prot;
        if (VirtualProtect(&ithunk->u1.Function, sizeof(void *),
                           PAGE_READWRITE, &prot)) {
          ithunk->u1.Function = (ULONGLONG) (ULONG_PTR) &elmer_win_c_setlocale;
          VirtualProtect(&ithunk->u1.Function, sizeof(void *), prot, &prot);
          patched++;
        }
      }
    }
  }
  return patched;
}

/* True if module base name matches "libgfortran*.dll" (case-insensitive). */
static int elmer_is_libgfortran(const char *name)
{
  size_t len;
  if (_strnicmp(name, "libgfortran", 11) != 0) return 0;
  len = strlen(name);
  return len >= 4 && _stricmp(name + len - 4, ".dll") == 0;
}

/* Enumerate loaded modules and patch every libgfortran*.dll. The module
 * enumeration APIs are resolved from kernel32 at runtime (they are the
 * K32*-prefixed forwarders present on Windows 7+), so no psapi link is
 * needed. Falls back to a few candidate sonames if enumeration is
 * unavailable. Safe no-op if libgfortran is not loaded. */
static void elmer_patch_libgfortran_setlocale(void)
{
  typedef BOOL  (WINAPI *EnumProcModules_fn)(HANDLE, HMODULE *, DWORD, LPDWORD);
  typedef DWORD (WINAPI *GetModuleBaseNameA_fn)(HANDLE, HMODULE, LPSTR, DWORD);

  HMODULE k32 = GetModuleHandleA("kernel32.dll");
  /* Cast through the generic function-pointer type to keep
   * -Wcast-function-type quiet about the FARPROC signature mismatch. */
  EnumProcModules_fn enum_modules = k32
      ? (EnumProcModules_fn) (void (*)(void)) GetProcAddress(k32, "K32EnumProcessModules") : NULL;
  GetModuleBaseNameA_fn base_name = k32
      ? (GetModuleBaseNameA_fn) (void (*)(void)) GetProcAddress(k32, "K32GetModuleBaseNameA") : NULL;

  if (enum_modules && base_name) {
    HMODULE mods[512];
    DWORD needed = 0;
    if (enum_modules(GetCurrentProcess(), mods, sizeof(mods), &needed)) {
      DWORD count = needed / sizeof(HMODULE);
      DWORD i;
      if (count > 512) count = 512;
      for (i = 0; i < count; i++) {
        char name[MAX_PATH];
        if (base_name(GetCurrentProcess(), mods[i], name, sizeof(name)) &&
            elmer_is_libgfortran(name)) {
          int n = elmer_patch_module_setlocale(mods[i]);
#ifdef ELMER_LOCALE_DEBUG
          fprintf(stderr, "[ElmerLocaleFix] %s: patched %d setlocale thunk(s)\n",
                  name, n);
#else
          (void) n;
#endif
        }
      }
      return;
    }
  }

  /* Fallback: probe the known/plausible sonames directly. */
  {
    static const char *const names[] = {
      "libgfortran-5.dll", "libgfortran-6.dll", "libgfortran-7.dll"
    };
    size_t i;
    for (i = 0; i < sizeof(names) / sizeof(names[0]); i++) {
      elmer_patch_module_setlocale(GetModuleHandleA(names[i]));
    }
  }
}
#endif /* _WIN32 */

/* Public entry, called from Fortran (ElmerSolver). Idempotent. */
void elmer_fix_numeric_locale(void)
{
  static int done = 0;
  if (done) return;
  done = 1;

  /* Elmer relies on "C" numeric formatting everywhere; make it explicit. */
  setlocale(LC_NUMERIC, "C");

#ifdef _WIN32
  elmer_patch_libgfortran_setlocale();
#endif
}
