#!/usr/bin/env python

# This is an automatically generated script to run your query
# to use it you will require the intermine python client.
# To install the client, run the following command from a terminal:
#
#     sudo easy_install intermine
#
# For further documentation you can visit:
#     http://intermine.readthedocs.org/en/latest/web-services/

# The line below will be needed if you are running this script with python 2.
# Python 3 will ignore it.
from __future__ import print_function

# The following two lines will be needed in every python script:
from intermine.webservice import Service
service = Service("https://yeastmine.yeastgenome.org/yeastmine/service")

# Get a new query on the class (table) you will be querying:
query = service.new_query("Gene")

# The view specifies the output columns
query.add_view(
    "phenotypeSummary", "phenotypes.experimentType", "phenotypes.mutantType",
    "phenotypes.observable", "phenotypes.strainBackground",
    "phenotypes.chemical", "phenotypes.condition", "phenotypes.details",
    "phenotypes.reporter", "name", "secondaryIdentifier", "briefDescription",
    "length", "symbol"
)

# Uncomment and edit the line below (the default) to select a custom sort order:
# query.add_sort_order("Gene.phenotypeSummary", "ASC")

for row in query.rows():
    print(row["phenotypeSummary"], row["phenotypes.experimentType"], row["phenotypes.mutantType"], \
        row["phenotypes.observable"], row["phenotypes.strainBackground"], \
        row["phenotypes.chemical"], row["phenotypes.condition"], row["phenotypes.details"], \
        row["phenotypes.reporter"], row["name"], row["secondaryIdentifier"], \
        row["briefDescription"], row["length"], row["symbol"])
 

