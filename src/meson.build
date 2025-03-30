# Compiler and build options
common_flags = []

# Detect the compiler and set flags accordingly
compiler = meson.get_compiler('cpp')

if compiler.get_id() == 'gcc' or compiler.get_id() == 'clang'
  extra_defines = ['-DNO_OLDNAMES']
  openmp_flag = get_option('openmp') ? ['-fopenmp'] : ['-DNO_OPENMP']
  zlib_link_flag = get_option('zlib') ? ['-lz'] : ['-DNO_ZLIB']
  flag_prefix = '-'
elif compiler.get_id() == 'msvc'
  extra_defines = ['/DNO_OLDNAMES']
  openmp_flag = get_option('openmp') ? ['/openmp'] : ['/DNO_OPENMP']
  zlib_link_flag = get_option('zlib') ? ['zlib.lib'] : ['/DNO_ZLIB']
  flag_prefix = '/'
else
  error('Unsupported compiler: ' + compiler.get_id())
endif

# Handle the MAX_SEQ option
max_seq_value = get_option('MAX_SEQ')
if max_seq_value > 0
  extra_defines += [flag_prefix + 'DMAX_SEQ=' + max_seq_value.to_string()]
endif

cc_flags = common_flags + openmp_flag + extra_defines
add_project_arguments(cc_flags, language: 'cpp')

deps = []
if get_option('zlib')
  # Fetch zlib dependency
  deps += [dependency('zlib', required: get_option('zlib'), fallback: ['zlib', 'zlib_dep'])]
endif


# Linker flags
link_static_flag = get_option('static') ? [flag_prefix + 'static'] : []
ld_flags = zlib_link_flag + link_static_flag

src_dir = join_paths(meson.project_source_root(), 'src')
message('Using source directory: ' + src_dir)

# Common source files
common_sources = [
  join_paths(src_dir, 'lib', 'aligners.cpp'),
  join_paths(src_dir, 'lib', 'cdhit-common.cpp'),
  join_paths(src_dir, 'lib', 'cdhit-utility.cpp'),
  join_paths(src_dir, 'lib', 'Options.ixx'),
  join_paths(src_dir, 'lib', 'ScoreMatrix.ixx'),
  join_paths(src_dir, 'lib', 'SequanceDB.ixx'),
  join_paths(src_dir, 'lib', 'WordTable.ixx'),
  join_paths(src_dir, 'lib', 'WorkingBuffer.ixx'),
]

# Source files for executables
sources = {
  'cd-hit': join_paths(src_dir, 'cdhit.cpp'),
  'cd-hit-2d': join_paths(src_dir, 'cdhit-2d.cpp'),
  'cd-hit-est': join_paths(src_dir, 'cdhit-est.cpp'),
  'cd-hit-est-2d': join_paths(src_dir, 'cdhit-est-2d.cpp'),
  'cd-hit-div': join_paths(src_dir, 'cdhit-div.cpp'),
  'cd-hit-454': join_paths(src_dir, 'cdhit-454.cpp'),
}

# Build targets
executables = []
foreach prog, src : sources
  exe = executable(
    prog,
    files(common_sources, src),
    install: true,
    link_args: ld_flags,
    install_dir: meson.project_source_root(),
    dependencies: deps,
  )
  executables += exe
endforeach
