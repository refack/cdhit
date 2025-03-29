# Compiler and build options
common_flags = []
openmp_flag = get_option('openmp') ? [] : ['-DNO_OPENMP']
zlib_flag = get_option('zlib') ? ['-DWITH_ZLIB'] : []

extra_defines = ['-DNO_OLDNAMES']
# Handle the MAX_SEQ option
max_seq_value = get_option('MAX_SEQ')
if max_seq_value > 0
  extra_defines += ['-DMAX_SEQ=' + max_seq_value.to_string()]
endif

cc_flags = common_flags + openmp_flag + zlib_flag + extra_defines
add_project_arguments(cc_flags, language: 'cpp')

# Linker flags
openmp_link_flag = get_option('openmp') ? [] : []
zlib_link_flag = get_option('zlib') ? ['-lz'] : []
link_static_flag = get_option('static') ? ['-static'] : []
ld_flags = openmp_link_flag + zlib_link_flag + link_static_flag

module_common = static_library(
  'module_common',
  'lib/Options.ixx',
  build_by_default : true,
)

common_sources = [
  'lib/cdhit-common.c++',
  'lib/cdhit-utility.c++',
  'lib/ScoreMatrix.cpp',
  'lib/aligners.cpp',
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
  exe = executable(prog, src, install: true, link_args: ld_flags, install_dir: meson.project_source_root(), link_with : module_common,)
  executables += exe
endforeach
