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
 * restore dance is redundant here. On Windows we therefore repoint
 * libgfortran-5.dll's imported setlocale to a no-op that reports "C": numeric
 * formatting stays correct (it already used "C") and the unsafe global locale
 * churn disappears. Only libgfortran's import is patched; Elmer's own and
 * MATC's setlocale() calls use a separate import table and are untouched.
 *
 * Called once from the ElmerSolver library entry point, on the main thread,
 * before any OpenMP region -> the plain static guard needs no locking.
 */

#include <locale.h>

#ifdef _WIN32
#include <windows.h>
#include <string.h>

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

/* Walk libgfortran-5.dll's import table and repoint every "setlocale" thunk
 * to elmer_win_c_setlocale. Safe no-op if the module or import is absent. */
static void elmer_patch_libgfortran_setlocale(void)
{
  HMODULE mod = GetModuleHandleA("libgfortran-5.dll");
  unsigned char *base;
  IMAGE_DOS_HEADER *dos;
  IMAGE_NT_HEADERS *nt;
  IMAGE_DATA_DIRECTORY dir;
  IMAGE_IMPORT_DESCRIPTOR *desc;

  if (!mod) return;
  base = (unsigned char *) mod;
  dos = (IMAGE_DOS_HEADER *) base;
  if (dos->e_magic != IMAGE_DOS_SIGNATURE) return;
  nt = (IMAGE_NT_HEADERS *) (base + dos->e_lfanew);
  if (nt->Signature != IMAGE_NT_SIGNATURE) return;
  dir = nt->OptionalHeader.DataDirectory[IMAGE_DIRECTORY_ENTRY_IMPORT];
  if (!dir.VirtualAddress) return;

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
        }
      }
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
