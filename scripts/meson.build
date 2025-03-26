# Install Perl scripts
perl_scripts = files(
    'FET.pl',
    'cd-hit-2d-para.pl',
    'cd-hit-clstr_2_blm8.pl',
    'cd-hit-div.pl',
    'cd-hit-para.pl',
    'clstr2tree.pl',
    'clstr2txt.pl',
    'clstr2xml.pl',
    'clstr_cut.pl',
    'clstr_list.pl',
    'clstr_list_sort.pl',
    'clstr_merge.pl',
    'clstr_merge_noorder.pl',
    'clstr_quality_eval.pl',
    'clstr_quality_eval_by_link.pl',
    'clstr_reduce.pl',
    'clstr_renumber.pl',
    'clstr_rep.pl',
    'clstr_reps_faa_rev.pl',
    'clstr_rev.pl',
    'clstr_select.pl',
    'clstr_select_rep.pl',
    'clstr_size_histogram.pl',
    'clstr_size_stat.pl',
    'clstr_sort_by.pl',
    'clstr_sort_prot_by.pl',
    'clstr_sql_tbl.pl',
    'clstr_sql_tbl_sort.pl',
    'make_multi_seq.pl',
    'plot_2d.pl',
    'plot_len1.pl',
)

python_exe = find_program('python3', 'python')
c = run_command(python_exe, 'glob_pl.py', check: true)
pl_validate = files(c.stdout().strip().split('\n'))

new_scripts = []
foreach f : pl_validate
    if not perl_scripts.contains(f)
        message('New Perl script found: ' + f.full_path())
        new_scripts += f
    endif
endforeach

if new_scripts.length() != 0
    msg = 'New Perl scripts found not explicitly mentioned:'
    foreach f : new_scripts
        msg += '\n  ' + f.full_path()
    endforeach
    error(msg)
endif

install_data(perl_scripts, install_dir: get_option('prefix') / 'bin')