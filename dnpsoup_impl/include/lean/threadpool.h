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


#ifndef LEAN_THREADPOOL_H
#define LEAN_THREADPOOL_H

#include <future>       // packaged_task
#include <functional>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <cstdint>
#include <queue>
#include <vector>
#include <memory>
#include <iostream>
#include <fstream>



namespace lean {
  inline namespace v1 {
    template<typename T>
    class ThreadPool {
    public:
        ThreadPool();
        ThreadPool(std::uint64_t n_workers);
        ~ThreadPool();

        ThreadPool<T>& add_task(std::function<T()>&&);

        ThreadPool<T>& run();
        std::vector<T> get_results();
    private:
        mutable std::mutex m_mu;
        mutable std::condition_variable m_cv;
        bool m_done;
        std::uint64_t m_worker_cnt;

        std::queue<std::function<T()>> m_job_q;
        std::vector<std::thread> m_workers;
        std::vector<std::future<std::vector<T>>> m_result_getters;
    };  // class ThreadPool
  } // v1
} // lean

#include "threadpool_impl.hpp"

#endif
