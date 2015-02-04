import logging
import os
import string


logger = logging.getLogger('report')


def create_rmd(summary_fn):
    root_path, fn = os.path.split(summary_fn)
    out_file = summary_fn.replace(".csv", "_re.csv")
    with open(summary_fn) as in_handle:
        with open(out_file, 'w') as out_handle:
            for line in in_handle:
                cols = line.strip().split(",")
                fix_line = ",".join([os.path.relpath(c, root_path) if os.path.exists(c) else c for c in cols])
                print >>out_handle, fix_line
    report_file = modify_report(root_path, fn)
    return out_file, report_file


def modify_report(summary_path, summary_fn):
    summary_path = os.path.abspath(summary_path)
    template = os.path.normpath(os.path.join(os.path.dirname(os.path.realpath(__file__)), "templates/report.rmd"))
    content = open(template).read()
    out_content = string.Template(content).safe_substitute({'path_abs': summary_path,
                                                            'path_summary': os.path.join(summary_path, summary_fn)})
    out_file = "ready_report.rmd"
    with open(out_file, 'w') as out_handle:
        print >>out_handle, out_content

    return out_file
