project('cdhit', 'cpp')

# Compiler and build options
cc = meson.get_compiler('cpp')
common_flags = ['-DNO_OLDNAMES']
openmp_flag = get_option('openmp') ? ['-fopenmp'] : ['-DNO_OPENMP']
zlib_flag = get_option('zlib') ? ['-DWITH_ZLIB'] : []
debug_flag = get_option('debug') ? ['-ggdb'] : ['-O2']

cc_flags = common_flags + openmp_flag + zlib_flag + debug_flag
add_project_arguments(cc_flags, language: 'cpp')

# Linker flags
openmp_link_flag = get_option('openmp') ? ['-fopenmp'] : []
zlib_link_flag = get_option('zlib') ? ['-lz'] : []
ld_flags = openmp_link_flag + zlib_link_flag

common_sources = [
  'cdhit-common.c++',
  'cdhit-utility.c++',
]

# Source files
sources = {
  'cd-hit': files(common_sources, 'cdhit.c++'),
  'cd-hit-2d': files(common_sources, 'cdhit-2d.c++'),
  'cd-hit-est': files(common_sources, 'cdhit-est.c++'),
  'cd-hit-est-2d': files(common_sources, 'cdhit-est-2d.c++'),
  'cd-hit-div': files(common_sources, 'cdhit-div.c++'),
  'cd-hit-454': files(common_sources, 'cdhit-454.c++'),
}

# Build targets
executables = []
foreach prog, src : sources
  exe = executable(prog, src, install: true, cpp_args: cc_flags, link_args: ld_flags)
  executables += exe
endforeach

subdir('scripts')
