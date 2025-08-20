# Compiler and build options
common_flags = []
extra_defines = ['-DNO_OLDNAMES']
ld_flags = []

# Detect the compiler and set flags accordingly
compiler = meson.get_compiler('cpp')

# Add AVX2 support
if get_option('enable_avx2')
  if compiler.get_id() == 'gcc' or compiler.get_id() == 'clang'
    common_flags += ['-mavx2']
  elif compiler.get_id() == 'msvc'
    common_flags += ['/arch:AVX2']
  endif
  extra_defines += ['-DUSE_AVX2']
endif

if compiler.get_id() == 'gcc' or compiler.get_id() == 'clang'
  if get_option('openmp')
    extra_defines += ['-DOPENMP', '-fopenmp']
  endif
  if get_option('zlib')
    ld_flags += ['-lz']
  else
    extra_defines += ['-DNO_ZLIB']
  endif
elif compiler.get_id() == 'msvc'
  extra_defines += [
#    '/WL',                   # Enable one-line diagnostics
    '/favor:INTEL64',        # Some special MSVC flags
  ]
  # extra_defines += ['/openmp']
  if meson.backend() == 'ninja'
    extra_defines += ['/reference', '"std=std.ifc"'] # Enable one-line diagnostics
  endif
else
  error('Unsupported compiler: ' + compiler.get_id())
endif

# Handle the MAX_SEQ option
max_seq_value = get_option('MAX_SEQ')
if max_seq_value > 0
  extra_defines += ['-D_MAX_SEQ_=' + max_seq_value.to_string()]
endif
extra_defines += ['-D_MAX_UAA_=4']

if compiler.get_id() == 'clang'
  cpp_flags = common_flags + extra_defines + ['-fmodules', '-fexperimental-library']
else
  cpp_flags = common_flags + extra_defines + ['-fmodules-ts']
endif
add_project_arguments(cpp_flags, language: 'cpp')
add_project_link_arguments(ld_flags, language: 'cpp')

deps = [
  dependency('zlib', static: get_option('static'), required: get_option('zlib')),
]
#include_directories: include_directories('./vcpkg_installed/x64-windows/include')
add_project_dependencies(deps, language: 'cpp')

src_dir = join_paths(meson.project_source_root(), 'src')
message('Using source directory: ' + src_dir)

# Common source files
modules_sources = [
  #
  #  join_paths(src_dir, 'lib', 'std.ixx'),
  #
  join_paths(src_dir, 'lib', 'common.ixx'),
  join_paths(src_dir, 'lib', 'gzfstream.ixx'),
  join_paths(src_dir, 'lib', 'Options.ixx'),
  join_paths(src_dir, 'lib', 'ScoreMatrix.ixx'),
  join_paths(src_dir, 'lib', 'Sequance.ixx'),
  join_paths(src_dir, 'lib', 'SequanceDB.ixx'),
  join_paths(src_dir, 'lib', 'WordTable.ixx'),
  join_paths(src_dir, 'lib', 'WorkingBuffer.ixx'),
  join_paths(src_dir, 'lib', 'WorkingParam.ixx'),
]
if compiler.get_id() == 'msvc'
  if meson.backend() == 'ninja'
    modules_sources += [join_paths(src_dir, 'microsoft', 'STL', 'std.ixx')]
  endif
endif
incs = []
if compiler.get_id() == 'msvc'
  incs = include_directories('C:/bin/dev/Microsoft/VS2022Community/VC/Tools/MSVC/14.43.34808/modules')
endif


# Source files for executables
sources = {
  # 'cd-hit': join_paths(src_dir, 'cdhit.cpp'),
  # 'cd-hit-2d': join_paths(src_dir, 'cdhit-2d.cpp'),
  'cd-hit-est': join_paths(src_dir, 'cdhit-est.cpp'),
  # 'cd-hit-est-2d': join_paths(src_dir, 'cdhit-est-2d.cpp'),
  # 'cd-hit-div': join_paths(src_dir, 'cdhit-div.cpp'),
  # 'cd-hit-454': join_paths(src_dir, 'cdhit-454.cpp'),
}

# Build targets
executables = []
foreach prog, src : sources
  exe = executable(
    prog,
    files(modules_sources, src),
    include_directories: incs,
    install: true,
    install_dir: meson.project_source_root(),
  )

  executables += exe
endforeach