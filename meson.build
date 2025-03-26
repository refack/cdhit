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

# Source files
sources = {
  'cd-hit': files('cdhit-common.c++', 'cdhit-utility.c++', 'cdhit.c++'),
  'cd-hit-2d': files('cdhit-common.c++', 'cdhit-utility.c++', 'cdhit-2d.c++'),
  'cd-hit-est': files('cdhit-common.c++', 'cdhit-utility.c++', 'cdhit-est.c++'),
  'cd-hit-est-2d': files('cdhit-common.c++', 'cdhit-utility.c++', 'cdhit-est-2d.c++'),
  'cd-hit-div': files('cdhit-common.c++', 'cdhit-utility.c++', 'cdhit-div.c++'),
  'cd-hit-454': files('cdhit-common.c++', 'cdhit-utility.c++', 'cdhit-454.c++'),
}

# Build targets
executables = []
foreach prog, src : sources
  exe = executable(prog, src, install: true, cpp_args: cc_flags, link_args: ld_flags)
  executables += exe
endforeach


