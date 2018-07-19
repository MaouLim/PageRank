#ifndef _PAGE_RANK_SPARE_MATRIX_F_H_
#define _PAGE_RANK_SPARE_MATRIX_F_H_

#include <string>
#include <fstream>
#include <sstream>
#include <set>
#include <vector>

#include <boost/smart_ptr/shared_array.hpp>

#include <primitives.h>
#include <cassert>
#include <mutex>

namespace tools {

	template <typename _IndexTp, typename _ValueTp>
	bool construct_mat_file(
		const std::string&                 destinate, 
		const std::set<
			std::pair<_IndexTp, _IndexTp>, 
			pair_less<_IndexTp, _IndexTp>
		>&                                 graph, 
		const std::vector<_IndexTp>&       count_in_degree
	) {

		typedef _IndexTp index_type;
		typedef _ValueTp value_type;

		if (graph.empty() || count_in_degree.empty()) { return false; }

		const size_t index_type_size = sizeof index_type;
		const size_t value_type_size = sizeof value_type;

		std::ofstream stream(destinate, std::ios::binary | std::ios::out);
		if (!stream.is_open()) { return false; }

		const size_t meta_size  = sizeof (size_t) * 2;
		const size_t index_size = index_type_size * (count_in_degree.size() + 1);
		const size_t pair_size  = value_type_size + index_type_size;
		const size_t data_size  = pair_size * graph.size();

		const size_t index_start = meta_size;
		const size_t data_start  = index_start + index_size;
		const size_t data_finish = data_start + data_size;

		stream.write(reinterpret_cast<const char*>(&data_start) , sizeof size_t);
		stream.write(reinterpret_cast<const char*>(&data_finish), sizeof size_t);

		stream.write(reinterpret_cast<const char*>(&data_start), index_type_size);

		size_t     count = 0;
		index_type row   = 0;
		for (const auto& each : graph) {
			const size_t data_offset = data_start + pair_size * count;
			stream.seekp(data_offset);

			stream.write(reinterpret_cast<const char*>(&each.second), index_type_size);
			value_type val = 1.0 / count_in_degree[each.second];
			stream.write(reinterpret_cast<const char*>(&val), value_type_size);

			if (row != each.first) {
				while (row < each.first) {
					++row;
					const size_t index_offset = index_start + index_type_size * row;
					stream.seekp(index_offset);

					static_assert(sizeof (size_t) <= index_type_size); // if index_type_size < size_t -> error!
					stream.write(reinterpret_cast<const char*>(&data_offset), index_type_size);
				}
			}
			++count;
		}

		size_t rest = index_start + index_type_size * (row + 1);
		while (rest < data_start) {
			stream.seekp(rest);
			stream.write(reinterpret_cast<const char*>(&data_finish), index_type_size);
			rest += index_type_size;
		}

		stream.close();
		return true;
	}

	/* 
	 * @note  Read from .txt, and convert it to a .data file which 
	 *        stores the spare matrix.
	 */
	template <typename _IndexTp, typename _ValueTp>
	bool convert(
		const std::string& source,
		const std::string& destinate
	) {
		typedef _IndexTp                      index_type;
		typedef std::pair<_IndexTp, _IndexTp> pair_type;
		typedef pair_less<_IndexTp, _IndexTp> less_op;

		std::ifstream stream(source, std::ios::in | std::ios::binary);
		if (!stream.is_open()) { return false; }

		size_t count = 0;
		bool count_completed = false;

		std::vector<index_type>      count_in_degree;
		std::set<pair_type, less_op> set;

		std::string line;
		while (std::getline(stream, line)) {
			if (!count_completed) {
				if (line.length() < 2) {
					count_in_degree.resize(count, static_cast<index_type>(0));
					count_completed = true;
					continue;
				}
				++count;
			}
			else {
				index_type row(0), col(0);
				std::stringstream sstream(line);
				sstream >> row >> col;
				if (row == col) { continue; }
				set.emplace(row, col);
				++count_in_degree[col];
			}
		}
		stream.close();

		return construct_mat_file<_IndexTp, _ValueTp>(destinate, set, count_in_degree);
	}

	template <
		typename _IndexTp, typename _ValTp
	>
	class spare_matrix_f {

	public:
		typedef _IndexTp index_type;
		typedef _ValTp   value_type;
		typedef _IndexTp dim_type;

		static const size_t index_type_size = sizeof index_type;
		static const size_t value_type_size = sizeof value_type;

	private:
		typedef spare_matrix_f<_IndexTp, _ValTp> self_type;
		typedef std::vector<index_type>          index_vec;
		typedef char                             byte_t;
		typedef boost::shared_array<byte_t>      buffer_type;

		static const size_t pair_size = index_type_size + value_type_size;

	public:
		explicit spare_matrix_f(const std::string& mat_file) : 
			m_stream(mat_file, std::ios::in | std::ios::binary)
		{
			assert(m_stream.is_open());
			_read_index_data();
		}

		~spare_matrix_f() { m_stream.close(); }

		dim_type dim() const { return m_index.size() - 1; }

		// todo not test enough
		value_type operator()(index_type row, index_type col) {
			buffer_type buffer;
			size_t buff_size = _read_row_data(row, buffer);

			buffer_iterator begin(buffer.get());
			buffer_iterator end(buffer.get() + buff_size);

			auto itr = std::lower_bound(
				begin, end, col, 
				[](const byte_t& first, index_type col) {
					return _column(&first) < col;
				}
			);

			size_t offset = std::distance(begin, itr);

			if (buff_size <= offset || col != _column(itr.base())) { return zero; }
			return _value(itr.base());
		}

		/* 
		 * @note  Apply a function [consumer] to each non-zero element 
		 *        of a row of matrix
		 */
		template <typename _Consumer>
		void apply_to_row(index_type row, _Consumer&& consumer) {
			buffer_type buffer;
			auto buff_size = _read_row_data(row, buffer);

			for (size_t i = 0; i < buff_size; i += pair_size) {
				const byte_t* ptr = &buffer[i];
				consumer(_column(ptr), _value(ptr));
			}
		}

	private:
		static value_type _value(const byte_t* data) {
			value_type val(0);
			memcpy(&val, data + index_type_size, value_type_size);
			return val;
		}

		static index_type _index(const byte_t* data) {
			index_type index(0);
			memcpy(&index, data, index_type_size);
			return index;
		}

		static index_type _column(const byte_t* data) {
			return _index(data);
		}

		size_t _read_row_data(index_type row, buffer_type& buffer) {
			const size_t bytes = m_index[row + 1] - m_index[row];
			if (0 == bytes) { return 0; }

			buffer.reset(new byte_t[bytes]);

			std::lock_guard<std::mutex> locker(m_mutex);

			m_stream.seekg(m_index[row]);
			m_stream.read(buffer.get(), bytes);

			return bytes;
		}

		void _read_index_data() {
			size_t index_start  = sizeof (size_t) * 2;
			size_t index_finish = 0;

			m_stream.seekg(0);
			m_stream.read(reinterpret_cast<byte_t*>(&index_finish), sizeof size_t);

			assert(
				index_start < index_finish && 0 == (index_finish - index_start) % index_type_size
			);

			const size_t buff_size = index_finish - index_start;
			buffer_type buffer(new byte_t[buff_size]);

			m_stream.seekg(index_start);
			m_stream.read(buffer.get(), buff_size);

			m_index.reserve(buff_size / index_type_size);
			for (size_t i = 0; i < buff_size; i += index_type_size) {
				m_index.emplace_back(_index(&buffer[i]));
			}
		}

		/* 
		 * @note  Inner class in order to adapt STL algorithm(std::lowwer_bound).
		 *        buffer_iterator is to iterate buffer_type objects, and it is
		 *        designed as a random asscess iterator.
		 */
		class buffer_iterator {

			typedef buffer_iterator self_type;
			typedef byte_t*         base_ptr;

		public:
			/* type traits */
			typedef byte_t                          value_type;
			typedef byte_t*                         pointer;
			typedef byte_t&                         reference;
			typedef std::random_access_iterator_tag iterator_category;
			typedef size_t                          difference_type;

			explicit buffer_iterator(base_ptr p) : m_ptr(p) { }

			base_ptr base() const { return m_ptr; }
			
			reference operator*() const { return *m_ptr; }
			pointer operator->() const { return &operator*(); }
			
			self_type& operator++() {
				m_ptr += pair_size;
				return *this;
			}
			
			self_type& operator++(int) {
				self_type tmp(m_ptr);
			  	m_ptr += pair_size;
			  	return *tmp;			}

			self_type& operator--() {
				m_ptr -= pair_size;
				return *this;
			}

			self_type& operator--(int) {
				self_type tmp(m_ptr);
				m_ptr -= pair_size;
				return *tmp;			}

			/* random access iterator constraints */
			reference operator[](difference_type n) const { return m_ptr[pair_size * n]; }

			self_type& operator+=(difference_type n) { m_ptr += pair_size * n; return *this; }
			self_type operator+(difference_type n) const { return self_type(m_ptr + pair_size * n); }

			self_type& operator-=(difference_type n) { m_ptr -= pair_size * n; return *this; }
			self_type operator-(difference_type n) const { return self_type(m_ptr - pair_size * n); }

			/* adapt to std::distance */
			difference_type operator-(const self_type& other) const { return m_ptr - other.m_ptr; }
			bool operator<(const self_type& other) const { return m_ptr < other.m_ptr; }
  
		private:
			base_ptr  m_ptr;
		};
			
	private:
		static const value_type zero;
		static const value_type one;

		mutable std::mutex m_mutex; /* mutex for file stream */
		std::ifstream      m_stream;
		index_vec          m_index;
	};

	template <typename _IndexTp, typename _ValTp>
	const _ValTp spare_matrix_f<_IndexTp, _ValTp>::zero(static_cast<_ValTp>(0));

	template <typename _IndexTp, typename _ValTp>
	const _ValTp spare_matrix_f<_IndexTp, _ValTp>::one(static_cast<_ValTp>(1));

	typedef spare_matrix_f<size_t, double> smatrixf;
}

#endif