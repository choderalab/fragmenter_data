import fragmenter

oemols = fragmenter.chemi.file_to_oemols('kinase_inhibitors.smi')
fragmenter.chemi.to_pdf(oemols, 'kinase_inhibitors.pdf')