# make sure test_filter.py is on PYTHONPATH so it can be found for import
from Filters import d3d_review_filter_1

import dream3d.simpl as simpl

dca = simpl.DataContainerArray()

# set up any other filters
print(f'Instantiating Filter D3DReviewTestFilter...')
py_filter = d3d_review_filter_1.D3DReviewTestFilter()
error_code, message = py_filter.execute(dca)

print('Done.')

