__author__ = 'ctroupin'

import pysocibclient

time_init = '2013-01-01T000000'
time_end = '2015-06-01T000000'

# Generate lists of platforms
socib_api = pysocibclient.ApiClient()
gliderlist = socib_api.list_platforms(instrument_type="glider")
gliderlist2 = socib_api.list_deployments(instrument_type="glider")

for glider in gliderlist:

    glider_opendap = glider.product_list[-1].opendap
    print glider_opendap


