#ifndef _PAGE_RANK_PRIMITIVES_H_
#define _PAGE_RANK_PRIMITIVES_H_

#include <utility>

namespace tools {
	
	template <typename _First, typename _Second>
	struct pair_less {

		typedef std::pair<_First, _Second> pair_type;

		bool operator()(const pair_type& a, const pair_type& b) const {
			if (a.first < b.first) { return true; }
			if (a.first == b.first && a.second < b.second) { return true; }
			return false;
		}
	};

}



#endif
