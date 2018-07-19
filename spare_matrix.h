/*
 * Created by Maou Lim on 2018-07-10
 */
#ifndef _PAGE_RANK_SPARE_MATRIX_H_
#define _PAGE_RANK_SPARE_MATRIX_H_

#include <cassert>
#include <sstream>
#include <algorithm>
#include <type_traits>

#include <vector.h>
#include <primitives.h>

namespace tools {

	struct high_precise { };
	struct mid_precise  { };
	struct low_precise  { };

	template <typename _ValTp>
	struct precise_trait {
		typedef low_precise precise_type;

		static const size_t precise_length = 1;
	};

	template <>
	struct precise_trait<float> {
		typedef mid_precise precise_type;
		typedef float       error_type;

		static const size_t precise_length = 6;
	};

	template <>
	struct precise_trait<double> {
		typedef high_precise precise_type;

		static const size_t precise_length = 16;
	};

	

	template <
		typename _IndexTp, 
		typename _ValTp,
		template <typename>
			class _VecTmpl = sequence
	>
	class spare_matrix {

	public:
		typedef _IndexTp index_type;
		typedef _ValTp   value_type;
		typedef _IndexTp dim_type;

		typedef _ValTp&       reference;
		typedef const _ValTp& const_reference;

		typedef _VecTmpl<value_type> vec_type;

	private:
		typedef spare_matrix<_IndexTp, _ValTp, _VecTmpl> self_type;

		typedef _VecTmpl<index_type> index_vec;
		typedef _VecTmpl<value_type> value_vec;

	public:
		enum class mat_type { ZERO, IDENTITY };
		enum class reshape_option { ZERO, IDENTITY, ADAPT };

		template <typename _ITp, typename _VTp, template <typename> class _VT, typename _VecTp>
		friend _VecTp operator*(const spare_matrix<_ITp, _VTp, _VT>&, const _VecTp&);

		template <typename _ITp, typename _VTp, template <typename> class _VT, typename _OthrVTp>
		friend spare_matrix<_ITp, _VTp, _VT> operator*(
			const spare_matrix<_ITp, _VTp, _VT>&, const _OthrVTp&
		);

		template <typename _ITp, typename _VTp, template <typename> class _VT, typename _OthrVTp>
		friend spare_matrix<_ITp, _VTp, _VT> operator*(
			const _OthrVTp&, const spare_matrix<_ITp, _VTp, _VT>&
		);

		spare_matrix(dim_type dim, mat_type type) : m_dim(dim) {
			if (mat_type::ZERO == type) {
				m_row_index.resize(m_dim + 1, index_type(0));
				return;
			}

			m_row_index.reserve(dim + 1);
			m_col_index.reserve(dim);
			m_values.reserve(dim);

			for (index_type i = 0; i < m_dim; ++i) {
				m_row_index.push_back(i);
				m_col_index.push_back(i);
				m_values.push_back(one);
			}

			m_row_index.push_back(dim);
		}

		spare_matrix(const self_type& other) :
			m_dim(other.m_dim),
			m_row_index(other.m_row_index),
			m_col_index(other.m_col_index),
			m_values(other.m_values) { }

		spare_matrix(self_type&& other) noexcept :
			m_dim(std::move(other.m_dim)),
			m_row_index(std::move(other.m_row_index)),
			m_col_index(std::move(other.m_col_index)),
			m_values(std::move(other.m_values)) { }

		self_type& operator=(const self_type& other) {  // NOLINT(misc-unconventional-assign-operator)
			if (this == &other) { return *this; }
				
			m_dim = other.m_dim;
			m_row_index = other.m_row_index;
			m_col_index = other.m_col_index;
			m_values = other.m_values;

			return *this;
		}

		self_type& operator=(self_type&& other) noexcept {
			if (this == &other) { return *this; }
				
			m_dim = std::move(other.m_dim);
			m_row_index = std::move(other.m_row_index);
			m_col_index = std::move(other.m_col_index);
			m_values = std::move(other.m_values);

			return *this;
		}

		const dim_type& dim() const { return m_dim; }

		void reshape(dim_type dim, reshape_option option = reshape_option::ADAPT) {
			switch (option) {
				case reshape_option::ZERO : {
					this->operator=(self_type(dim, mat_type::ZERO));
					break;
				}
				case reshape_option::IDENTITY : {
					this->operator=(self_type(dim, mat_type::IDENTITY));
					break;
				}
				case reshape_option::ADAPT : {
					if (dim == m_dim) { return; }

					dim_type min = std::min(dim, m_dim);
					self_type tmp(dim, mat_type::ZERO);

					// todo  this can be optimized
					for (index_type row(0); row < min; ++row) {
						for (index_type col(0); col < min; ++col) {
							tmp.set(row, col, this->operator()(row, col));
						}
					}

					this->operator=(std::move(tmp));
					return;
				}
				default: { assert(false); }
			}
		}

		const_reference operator()(index_type row, index_type col) const {

			auto itr = std::lower_bound(
				m_col_index.begin() + m_row_index[row], 
				m_col_index.begin() + m_row_index[row + 1], 
				col
			);

			index_type offset = 
				static_cast<index_type>(std::distance(m_col_index.begin(), itr));

			if (m_row_index[row + 1] <= offset || col != *itr) { return zero; }
			
			return m_values[offset];
		}

		template <typename _OtherValTp>
		void set(index_type row, index_type col, const _OtherValTp& val) {
			auto itr = std::lower_bound(
				m_col_index.begin() + m_row_index[row],
				m_col_index.begin() + m_row_index[row + 1],
				col
			);

			index_type offset =
				static_cast<index_type>(std::distance(m_col_index.begin(), itr));

			if (m_row_index[row + 1] <= offset || col != *itr) { 
				/* not found -> mat(row, col) == zero. */
				if (eps <= std::abs(val)) {
					for (index_type i(row + 1); i < m_dim + 1; ++i) {
						++m_row_index[i];
					}
					m_col_index.insert(m_col_index.begin() + offset, col);
					m_values.insert(m_values.begin() + offset, static_cast<value_type>(val));
				}
			}
			/* found -> mat(row, col) != zero */
			else if (std::abs(val) < eps) {
				/* set mat(row, col) zero */
				for (index_type i(row + 1); i < m_dim + 1; ++i) {
					--m_row_index[i];
				}
				m_col_index.erase(m_col_index.begin() + offset);
				m_values.erase(m_values.begin() + offset);
			}
			else { m_values[offset] = val; }
		}

		/* 
		 * @note  ignore zeros 
		 */
		template <typename _BinaryOp>
		value_type accumulation_of_row(
				index_type row, _BinaryOp&& op
		) const {
			value_type sum(static_cast<value_type>(0));

			for (index_type i(m_row_index[row]); i != m_row_index[row + 1]; ++i) {
				sum = op(sum, m_values[i]);
			}

			return sum;
		}

		/*
		 * @note  ignore zeros
		 */
		template <typename _BinaryOp>
		value_type accumulation_of_col(
				index_type col, _BinaryOp&& op
		) const {
			value_type sum(static_cast<value_type>(0));

			auto cur = m_col_index.begin();

			while (m_col_index.end() != cur) {
				cur = std::find(cur, m_col_index.end(), col);
				if (m_col_index.end() == cur) { break; }
				index_type offset = std::distance(m_col_index.begin(), cur++);
				sum = op(sum, m_values[offset]);
			}
			
			return sum;
		}

		template <typename _Consumer>
		void apply_to_row(index_type row, _Consumer&& consumer) {
			for (index_type i(m_row_index[row]); i != m_row_index[row + 1]; ++i) {
				m_values[i] = consumer(m_col_index[i], m_values[i]);
			}
		}

		template <typename _Consumer>
		void apply_to_col(index_type col, _Consumer&& consumer) {
			auto cur = m_col_index.begin();

			while (m_col_index.end() != cur) {
				cur = std::find(cur, m_col_index.end(), col);
				if (m_col_index.end() == cur) { break; }
				index_type offset = std::distance(m_col_index.begin(), cur++);
				m_values[offset] = consumer(m_values[offset]);
			}
		}

		vec_type row_vec(index_type row) const {
			vec_type row_vec(m_dim, zero);

			index_type i(m_row_index[row]);
			while (i < m_row_index[row + 1]) {
				row_vec[m_col_index[i]] = m_values[i];
				++i;
			}

			return row_vec;
		}

		std::string to_smat_str() const {
			std::stringstream stream;
			stream.precision(precise_trait<value_type>::precise_length);

			stream << m_dim << std::endl;

			for (auto& each : m_row_index) {
				stream << each << ' ';
			}
			stream << std::endl;

			for (auto& each : m_col_index) {
				stream << each << ' ';
			}
			stream << std::endl;

			for (auto& each : m_values) {
				stream << each << ' ';
			}

			return stream.str();
		}

		std::string to_str() const {
			std::stringstream stream;
			stream.precision(precise_trait<value_type>::precise_length + 1);

			for (index_type i = 0; i < m_dim; ++i) {
				for (index_type j = 0; j < m_dim; ++j) {
					stream << this->operator()(i, j) << '\t';
				}
				stream << std::endl;
			}

			return stream.str();
		}

	public:
		static const value_type zero;
		static const value_type eps;
		static const value_type one;

	private:
		dim_type  m_dim;

		index_vec m_row_index;
		index_vec m_col_index;

		value_vec m_values;
	};

	template <typename _IndexTp, typename _ValTp, template <typename> class _VecTmpl>
	const _ValTp spare_matrix<_IndexTp, _ValTp, _VecTmpl>::zero(static_cast<_ValTp>(0));

	template <typename _IndexTp, typename _ValTp, template <typename> class _VecTmpl>
	const _ValTp spare_matrix<_IndexTp, _ValTp, _VecTmpl>::eps(static_cast<_ValTp>(1e-8));

	template <typename _IndexTp, typename _ValTp, template <typename> class _VecTmpl>
	const _ValTp spare_matrix<_IndexTp, _ValTp, _VecTmpl>::one(static_cast<_ValTp>(1));

	typedef spare_matrix<size_t, double, sequence> smatrix;

	template <
		typename _ITp, 
		typename _VTp, 
		template <typename> 
			class _VT, 
		typename _VecTp
	>
	_VecTp operator*(
		const spare_matrix<_ITp, _VTp, _VT>& mat, const _VecTp& vec
	) {
		typedef _VecTp                        vec_t;
		typedef spare_matrix<_ITp, _VTp, _VT> smat_t;
		typedef typename smat_t::index_type   idx_t;
		typedef _VTp                          val_t;
 
		//static_assert(std::is_same<val_t, typename vec_t::value_type>::value);
		assert(vec.size() == mat.m_dim);
 
		vec_t result(mat.m_dim, smat_t::zero);
 
		for (idx_t i = 0; i < mat.m_dim; ++i) {
			for (idx_t j = mat.m_row_index[i]; j < mat.m_row_index[i + 1]; ++j) {
				result[i] += mat.m_values[j] * vec[mat.m_col_index[j]];
			}
		}
 
		return result;
	}

	template <
		typename _ITp,
		typename _VTp,
		template <typename>
			class _VT,
		typename _OthrVTp
	>
	spare_matrix<_ITp, _VTp, _VT> operator*(
		 const spare_matrix<_ITp, _VTp, _VT>& mat, const _OthrVTp& scalar
	) {
		typedef spare_matrix<_ITp, _VTp, _VT> smat_t;
		typedef typename smat_t::index_type   idx_t;
		typedef _VTp                          val_t;

		if (val_t(0) == scalar) { return smat_t(mat.m_dim, smat_t::mat_type::ZERO); }

		smat_t result(mat.m_dim, smat_t::mat_type::ZERO);

		for (idx_t row(0); row < mat.m_dim; ++row) {
			for (
				idx_t i(mat.m_row_index[row]); 
				i < mat.m_row_index[row + 1];
				++i
			) {
				val_t val = mat.m_values[i] * scalar;
				result.set(row, mat.m_col_index[i], val);
			}
		}

		return result;
	}

	template <
		typename _ITp,
		typename _VTp,
		template <typename>
			class _VT,
		typename _OthrVTp
	>
	spare_matrix<_ITp, _VTp, _VT> operator*(
		const _OthrVTp& scalar, const spare_matrix<_ITp, _VTp, _VT>& mat
	) { return tools::operator*(mat, scalar); }


	template <typename _IndexTp, typename _ValTp>
	bool load_graph(
		const std::string&              source,
		spare_matrix<_IndexTp, _ValTp>& graph
	) {
		typedef spare_matrix<_IndexTp, _ValTp> mat_t;
		typedef std::set<
			std::pair<_IndexTp, _IndexTp>,
			pair_less<_IndexTp, _IndexTp>
		> set_type;

		std::ifstream stream(source, std::ios::in);
		if (!stream.is_open()) { return false; }

		size_t count = 0;
		bool count_completed = false;

		set_type            set;
		std::vector<size_t> count_in_degree;

		std::string line;
		while (std::getline(stream, line)) {
			if (!count_completed) {
				if (line.empty()) {
					count_completed = true;
					count_in_degree.resize(count, 0);
					continue;
				}
				++count;
			}
			else {
				size_t row = 0, col = 0;
				std::stringstream sstream(line);
				sstream >> row >> col;
				if (row == col) { continue; }
				set.emplace(row, col);
				++count_in_degree[col];
			}
		}
		stream.close();

		mat_t tmp(count, mat_t::mat_type::ZERO);
		for (auto& each : set) {
			tmp.set(each.first, each.second, mat_t::one / count_in_degree[each.second]);
		}

		graph = std::move(tmp);
		return true;
	}
}

namespace std {

	/* 
	 * @note  expansion for std::operator<<
	 */
	template <typename _IndexTp, typename _ValTp>
	std::ostream& operator<<(
		std::ostream& stream, const tools::spare_matrix<_IndexTp, _ValTp>& mat
	) {
		return stream << mat.to_str();
	}
}

#endif