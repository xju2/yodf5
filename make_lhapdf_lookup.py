#!/usr/bin/env python

import argparse
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="make a dictionary for LHAPDF")
    add_arg = parser.add_argument
    add_arg("file_name", help='input csv file download from LHAPDF page')
    add_arg('-o', '--output', help='output name', default='lhapdf_lookup.py')
    add_arg('--cxx-output', help='c++ output name', default='lhapdf_lookup.hpp')
    args = parser.parse_args()

    ihead = 0
    out = "_dictionary = {\n"
    cxx_out = "#ifndef __YODF5_LHAPDF_LOOKUP_H__\n#define __YODF5_LHAPDF_LOOKUP_H__\n\n"
    cxx_out += "#include <map>\n#include <string>\n\nconst std::map<int, std::string> LHAPDF_DICT {\n"

    allowed_pdf = [
        'CT14nnlo', 'MMHT2014nnlo68cl', 'NNPDF30_nnlo_as_0118', 'NNPDF30_nnlo_as_0117', 'NNPDF30_nnlo_as_0119',
        "PDF4LHC15_nnlo_30_pdfas"
        ]
    with open(args.file_name, 'r') as f:
        for line in f:
            if ihead == 0:
                ihead += 1
                continue
            ihead += 1
            items = [x.strip('"') for x in line[:-1].split(',')]
            pdf_id_base = int(items[0])
            n_variations = int(items[2])
            for ivar in range(n_variations):
                pdf_id = pdf_id_base + ivar
                pdf_name = items[1]
                if pdf_name in allowed_pdf:
                    out += '    {:d}: "{}",\n'.format(pdf_id, items[1])
                    cxx_out += '    {{{:d}, "{}"}},\n'.format(pdf_id, items[1])
    out += '}\n'
    out += """\n
def pset_name(pid):
    return _dictionary[pid]
    """
    cxx_out += "};\n\n"
    
    with open(args.output, 'w') as f:
        f.write(out)

    cxx_out += """
std::string pset_name(int pid) {
    return LHAPDF_DICT.at(pid);
}\n

#endif
    """
    with open(args.cxx_output, 'w') as f:
        f.write(cxx_out)