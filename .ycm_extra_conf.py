
def Settings( **kwargs ):
  return {
    'flags': ['-x', 'c++', '-Wall', '-Wextra', '-Werror',
              '-isystem',
              '-std=++17',
              '-Wc++17-extensions',
              '-I./dnpsoup_impl/include', '-I./matrix/matrix_impl/include',
              '-I./dnpsoup_impl',
              '-I./matrix/matrix_impl',
              '-I/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/Headers',
              '-I./build',
              '-I./build/googletest/googletest-src/googletest/include'],
  }
