sed -e '/\#if defined(UCFUNNAME)/{;N;N;N;N;s/#if defined(UCFUNNAME)\nEIOFC(\([^\)]*\))\n\#else\n[ ]*EIOFC(\([^\)]*\))\n\#endif/FC_NAME_(\2,\1)/g}' eio_api.cpp
