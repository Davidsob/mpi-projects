FILE(REMOVE_RECURSE
  "libmodels.pdb"
  "libmodels.a"
)

# Per-language clean rules from dependency scanning.
FOREACH(lang)
  INCLUDE(CMakeFiles/models.dir/cmake_clean_${lang}.cmake OPTIONAL)
ENDFOREACH(lang)
