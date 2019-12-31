
def Settings( **kwargs ):
  return {
    'flags': ['-x', 'c++', '-Wall', '-Wextra', '-Werror',
              '-std=++17',
              '-Wc++17-extensions',
              '-I./dnpsoup_impl/include', '-I./matrix/matrix_impl/include',
              '-I./dnpsoup_impl',
              '-I./build/googletest/googletest-src/googletest/include',
              '-I./build'
              ]
  }
