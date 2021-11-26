#!/bin/bash
#Files provided have uneven number of spaces making it difficult to distinguish columns
find data -name '*.tab' -exec sed -ie 's/  */;/g' {} \; #substitutes any-number-of-spaces with a semicolon
find data -name '*.tab' -exec sed -ie 's/^;//g' {} \;   #removes leading semicolons
find data -name '*.tab' -exec sed -ie 's!Energy\[;kJ/mol;\]!Energy\[kJ/mol\]!g' {} \; #fixes header - there are inconsistent spaces in column description
