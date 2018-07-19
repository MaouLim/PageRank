/*
* Created by Maou Lim on 2018-07-13
*/
#ifndef _PAGE_RANK_SOVLER_H_
#define _PAGE_RANK_SOVLER_H_

#include <future>

#include <spare_matrix.h>
#include <spare_matrix_f.h>

namespace page_rank {
	
	template <
		typename _Result, typename _Params
	>
	struct solver {
	
		typedef _Result result_type;
		typedef _Params params_type;

		solver() = default;
		virtual ~solver() = default;

		virtual result_type solve(params_type&&) { return result_type(); }
		virtual result_type solve(params_type&) { return result_type(); }
	};

	class page_rank_solver : 
		solver<std::pair<size_t ,tools::vector>, tools::smatrix> {
		
		typedef std::pair<size_t, tools::vector>  pair_type;
		typedef solver<pair_type, tools::smatrix> base_type;
		typedef page_rank_solver                  self_type;

	public:
		typedef typename base_type::result_type result_type;
		typedef typename base_type::params_type params_type;

		explicit page_rank_solver(
			double alpha, size_t max_itr = 100, double err = 1e-8
		) : m_alpha(alpha), m_max_itr(max_itr), m_err(err) { }

		result_type solve(params_type& matrix) override {

			const size_t dim = matrix.dim();

			for (size_t row = 0; row < dim; ++row) {
				matrix.apply_to_row(
					row, [this](size_t, double item) { return item * m_alpha; }
				);
			}

			auto bias = (1.0 - m_alpha) / dim;
			auto prev = tools::generate_rand_vec<tools::vector>(dim);
			auto cur  = tools::vector(dim);

			size_t count = 0;

			std::vector<std::future<void>> futures(dim);

			while (count < m_max_itr) {

				for (size_t row = 0; row < dim; ++row) {
					futures[row] = std::async(
						std::launch::async, 
						&_product_task, 
						std::ref(matrix), 
						bias, 
						row, 
						std::ref(cur), 
						std::ref(prev)
					);
				}
				for (auto& each : futures) { each.wait(); }

				tools::normalize(cur);
				if (tools::euler_distance(cur, prev) < m_err) { break; }

				cur.swap(prev);
				++count;
			}

			return std::make_pair(count, cur);
		}

		double alpha() const { return m_alpha; }
		void set_alpha(double alpha) { m_alpha = alpha; }

		size_t max_iterations() const { return m_max_itr; }
		void set_max_iterations(size_t max_itr) { m_max_itr = max_itr; }

		double error() const { return m_err; }
		void set_error(double err) { m_err = err; }

	private:
		static void _product_task(
			params_type&   matrix, 
			double         bias, 
			size_t         row, 
			tools::vector& cur,
			tools::vector& prev
		) {
			matrix.apply_to_row(
				row,
				[bias, row, &cur, &prev](size_t col, double item) {
					cur[row] += (item + bias) * prev[col];
					return item;
				}
			);
		}

	private:
		double m_alpha;
		size_t m_max_itr;
		double m_err;
	};

	class page_rank_solver_f : 
		solver<std::pair<size_t, tools::vector>, tools::smatrixf> {

		typedef std::pair<size_t, tools::vector>   pair_type;
		typedef solver<pair_type, tools::smatrixf> base_type;
		typedef page_rank_solver_f                 self_type;

	public:
		typedef typename base_type::result_type result_type;
		typedef typename base_type::params_type params_type;

		explicit page_rank_solver_f(
			double alpha, size_t max_itr = 100, double err = 1e-8
		) : m_alpha(alpha), m_max_itr(max_itr), m_err(err) { }

		/*
		 * @param matrix  Gm matrix.
		 * @ret           The iterations and the target vector.
		 */
		result_type solve(params_type& matrix) override {

			const size_t dim = matrix.dim();

			auto bias = (1.0 - m_alpha) / dim;
			auto prev = tools::generate_rand_vec<tools::vector>(dim);
			auto cur = tools::vector(dim);

			size_t count = 0;

			std::vector<std::future<void>> futures(dim);

			while (count < m_max_itr) {

				for (size_t row = 0; row < dim; ++row) {

					/* 
					 * @note  Start an async-task to calculate the product of one row 
					 *        of [matrix] with vector [prev], the result will be stored
					 *        in vector [cur]
					 */
					futures[row] = std::async(
						std::launch::async,
						&_product_task,
						std::ref(matrix),
						m_alpha,
						bias,
						row,
						std::ref(cur),
						std::ref(prev)
					);

				}
				for (auto& each : futures) { each.wait(); }

				tools::normalize(cur);
				if (tools::euler_distance(cur, prev) < m_err) { break; }

				cur.swap(prev);
				++count;
			}

			return std::make_pair(count, cur);
		}

		double alpha() const { return m_alpha; }
		void set_alpha(double alpha) { m_alpha = alpha; }

		size_t max_iterations() const { return m_max_itr; }
		void set_max_iterations(size_t max_itr) { m_max_itr = max_itr; }

		double error() const { return m_err; }
		void set_error(double err) { m_err = err; }

	private:
		static void _product_task(
			params_type&   matrix,
			double         alpha,
			double         bias,
			size_t         row,
			tools::vector& cur,
			tools::vector& prev
		) {
			size_t last_nonzero = 0;
			matrix.apply_to_row(
				row,
				[alpha, bias, row, &last_nonzero, &cur, &prev](size_t col, double item) {
					//assert(col < prev.size());
					cur[row] += (alpha * item + bias) * prev[col];
					if (0.0 == bias) { return; }
					for (auto i = last_nonzero + 1; i < col; ++i) {
						cur[row] += bias * prev[i];
					}
					last_nonzero = col;
				}
			);
		}

	private:
		double m_alpha;
		size_t m_max_itr;
		double m_err;
	};

}

#endif