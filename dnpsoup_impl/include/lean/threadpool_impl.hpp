/*BSD 3-Clause License

Copyright (c) 2019, Chen Yang, yangcnju@gmail.com
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its
   contributors may be used to endorse or promote products derived from
   this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

namespace lean {
	inline namespace v1 {
		template<typename T>
		ThreadPool<T>::ThreadPool()
			: m_done(false), m_worker_cnt(std::thread::hardware_concurrency())
		{
		}

		template<typename T>
		ThreadPool<T>::ThreadPool(std::uint64_t n_workers)
			: m_done(false), m_worker_cnt(n_workers)
		{
		}

		template<typename T>
		ThreadPool<T>::~ThreadPool()
		{
			{
				std::lock_guard<std::mutex> lk(m_mu);
				m_done = true;
				while (!m_job_q.empty()) {
					m_job_q.pop();
				}
			}
			
			for (auto& t : m_workers) {
				t.join();
			}
		}

		template<typename T>
		ThreadPool<T>& ThreadPool<T>::add_task(std::function<T()>&& task)
		{
			{
				std::lock_guard<std::mutex> lk(m_mu);
				m_job_q.push(std::move(task));
			}
			
			return *this;
		}

		template<typename T>
		ThreadPool<T>& ThreadPool<T>::run()
		{
			for (std::uint64_t i = 0; i < m_worker_cnt; ++i) {
				std::packaged_task<std::vector<T>()> worker([this]() {
					std::vector<T> results;
					while (true) {
						std::unique_lock<std::mutex> lk(m_mu);
						m_cv.wait(lk, [this]() {
							return this->m_done || (!m_job_q.empty());
						});

						if (this->m_done && m_job_q.empty()) {
							m_cv.notify_one();
							break;
						}

						std::function<T()> task = std::move(m_job_q.front());
						m_job_q.pop();
						m_cv.notify_one();
						lk.unlock();

						T res = task();
						results.push_back(std::move(res));
					}
					return results;
				});
				m_result_getters.emplace_back(worker.get_future());
				m_workers.emplace_back(std::move(worker));
			}
      m_cv.notify_one();

			return *this;
		}

		template<typename T>
		std::vector<T> ThreadPool<T>::get_results()
		{
			std::vector<T> results;
			{
				std::unique_lock<std::mutex> lk(m_mu);
				m_cv.wait(lk, [this]() {
					return m_job_q.empty();
					});
				m_done = true;
				m_cv.notify_one();
			}

			for (auto& result_getter : m_result_getters) {
				auto res = result_getter.get();
				results.insert(results.end(), res.begin(), res.end());
			}
			return results;
		}
	}
}
