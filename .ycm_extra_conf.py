flags = [
    '-Wall',
    '-Wextra',
    '-Werror',
    '-std=c++17',
    '-I./dnpsoup_impl/include',
    '-I./matrix/matrix_impl/include'
]

def Settings( **kwargs ):
  return {
    'flags': [ '-x', 'c++', '-Wall', '-Wextra', '-Werror',
              '-I./dnpsoup_impl/include', '-I./matrix/matrix_impl/include'],
  }
