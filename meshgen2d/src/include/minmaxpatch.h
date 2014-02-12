#ifndef MINGW32
#ifndef MINMAXPATCH_H
#define MINMAXPATCH_H

namespace std
{

#undef min
#undef max

	template <class _Tp>
	inline const _Tp& min(const _Tp& __a, const _Tp& __b) {
	  return __b < __a ? __b : __a;
	}

	template <class _Tp>
	inline const _Tp& max(const _Tp& __a, const _Tp& __b) {
	  return  __a < __b ? __b : __a;
	}

};

#endif
#endif
