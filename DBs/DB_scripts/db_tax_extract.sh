#!/bin/bash

# usage: db_tax_extract.sh 

cat $1 \
  | grep '^>' \
  | awk -F'|' '
      { id = substr($1,2);    # remove leading ">"
        m = split($4,a,",");
        printf("%s", id);
        for(i=1;i<=9;i++){
          if(i<=m) { val=a[i]; gsub(/^ +| +$/,"",val); if(val=="") val="NA" }
          else val="NA";
          printf("\t%s", val);
        }
        printf("\n");
      }' \
  > taxa_table.tsv

