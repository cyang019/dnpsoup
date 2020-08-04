#ifndef LEAN_TASKQUEUE_H
#define LEAN_TASKQUEUE_H

#include <mutex>
#include <thread>
#include <future>
#include <functional>
#include <queue>
#include <memory>
#include <iostream>


namespace lean {
  template<typename T>
  class TaskQueue {
  public:
    TaskQueue();
    ~TaskQueue();

    TaskQueue<T>& add_task(std::function<T()>&& t);
    std::vector<T> run(size_t n_workers = std::thread::hardware_concurrency());
  private:
    //mutable std::mutex mu_;
    std::queue<std::function<T()>> job_q_;
  };

  template<typename T>
  TaskQueue<T>::TaskQueue()
  {}

  template<typename T>
  TaskQueue<T>::~TaskQueue()
  {
    while(!job_q_.empty()) {
      job_q_.pop();
    }
  }

  template<typename T>
  TaskQueue<T>& TaskQueue<T>::add_task(std::function<T()>&& task)
  {
    job_q_.push(std::move(task));
    return *this;
  }

  template<typename T>
  std::vector<T> TaskQueue<T>::run(size_t n_workers)
  {
    size_t n_busy = 0;
    bool done = false;
    std::vector<std::future<T>> result_getters;
    std::vector<std::thread> threads;

    auto execTask = [](std::function<T()> &&task) {
      auto td_task = std::move(task);
      T res = td_task();
      return res;
    };

    while(n_busy < n_workers) {
      if(job_q_.empty()) {
        break;
      }
#ifndef NDEBUG
      std::cout << "job_q_.size(): " << job_q_.size() << std::endl;
#endif
      std::function<T()> task = std::move(job_q_.front());
      job_q_.pop();
      std::packaged_task<T(std::function<T()> &&)> pkg_task(execTask);
      result_getters.emplace_back(pkg_task.get_future());
      std::thread task_td(std::move(pkg_task), std::move(task));
      threads.push_back(std::move(task_td));
      ++n_busy;
    }
#ifndef NDEBUG
    std::cout << "tasks running..." << std::endl;
#endif
    std::vector<bool> threads_joined(threads.size(), false);
    while (true) {
      if(done) {
        if (job_q_.empty())
          break;
      }
#ifndef NDEBUG
      std::cout << "job_q_.size() in the 2nd while loop: " << job_q_.size() << std::endl;
#endif
      for(size_t i = 0; i < threads.size(); ++i) {
        if(threads_joined[i]) {
          continue;
        }
        else if(threads[i].joinable()) {
#ifndef NDEBUG
          std::cout << "trying to join one thread." << std::endl;
          std::cout << "n_busy: " << n_busy << std::endl;
#endif
          threads[i].join();
#ifndef NDEBUG
          std::cout << "one thread joined." << std::endl;
#endif
          threads_joined[i] = true;
          --n_busy;
#ifndef NDEBUG
          std::cout << "n_busy: " << n_busy << std::endl;
#endif
          /// if a worker finishes, check if there is more work to do
          if(n_busy < n_workers && !job_q_.empty()) { ///< if idle worker exists
#ifndef NDEBUG
            std::cout << "Trying to retrieve job..." << std::endl;
            std::cout << "job_q_.size(): " << job_q_.size() << std::endl;
#endif
            std::function<T()> task = std::move(job_q_.front());
            job_q_.pop();
            std::packaged_task<T(std::function<T()> &&)> pkg_task(execTask);
            result_getters.emplace_back(pkg_task.get_future());
            std::thread task_new(std::move(pkg_task), std::move(task));
#ifndef NDEBUG
            std::cout << "trying to replace one thread." << std::endl;
            std::cout << "n_busy: " << n_busy << std::endl;
#endif
            threads[i] = std::move(task_new);
            threads_joined[i] = false;
#ifndef NDEBUG
            std::cout << "one thread replaced..." << std::endl;
#endif
            ++n_busy;
#ifndef NDEBUG
            std::cout << "n_busy: " << n_busy << std::endl;
#endif
          }
        }
      }
      if(n_busy == 0 && job_q_.empty()) {
        done = true;
      }
#ifndef NDEBUG
      std::cout << "x" << std::flush;
#endif
      std::this_thread::sleep_for(std::chrono::milliseconds(100));
    }

    std::vector<T> results;
    for(auto& result_getter : result_getters) {
      auto res = result_getter.get();
      results.push_back(res);
    }
    return results;
  }
} // namespace lean


#endif
