dnl
dnl Allow the user disable support for command line editing using GNU
dnl readline.
dnl
dnl derived from OCTAVE_READLINE...
dnl 
AC_DEFUN([ACX_ENABLE_READLINE], [
USE_READLINE=true
LIBREADLINE=
AC_ARG_ENABLE(readline,
    [  --enable-readline       use readline library (default is yes)],
    [if test "$enableval" = no; then
       USE_READLINE=false
       warn_readline="command editing and history features require GNU Readline"
     fi])
  if $USE_READLINE; then
    AC_CHECK_LIB(readline, rl_set_keyboard_input_timeout, [
      LIBREADLINE="-lreadline"
      LIBS="$LIBREADLINE $LIBS"

       
      dnl figure out where the they have hidden the header...
      AC_CHECK_HEADERS([readline.h],[],
	[
	   AC_CHECK_HEADERS([readline/readline.h],[],USE_READLINE=false)
	])

      if test $USE_READLINE == true; then
        AC_DEFINE(USE_READLINE, 1, [Define to use the readline library.])
      else 
        AC_MSG_WARN([Headers not found, disabling readline support.])
      fi
    ], [
      AC_MSG_WARN([I need GNU Readline 4.2 or later... trying to do without...])
    ],"-lcurses")
  fi
  AC_SUBST(LIBREADLINE)
])


AC_DEFUN([ACX_GET_TERM], [
acx_found_termlib=no
for termlib in ncurses curses termcap terminfo termlib; do
  AC_CHECK_LIB(${termlib}, tputs, [TERMLIBS="${TERMLIBS} -l${termlib}"])
  case "${TERMLIBS}" in
    *-l${termlib}*)
      LIBS="$TERMLIBS $LIBS"
      AC_SUBST(TERMLIBS)
      acx_found_termlib=yes
      break
    ;;
  esac
done

if test "$acx_found_termlib" = no; then
  warn_termlibs="I couldn't find -ltermcap, -lterminfo, -lncurses, -lcurses, o\
r -ltermlib!"
  AC_MSG_WARN($warn_termlibs)
fi
])
