project('cdhit', 'cpp', version : '25.0.0', default_options: ['cpp_std=c++latest'])
add_project_arguments('-DCDHIT_VERSION="' + meson.project_version() + '"', language: 'cpp')
# add_project_arguments('-std=c++2c', language: 'cpp')

subdir('src')
subdir('scripts')