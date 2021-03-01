# tidyrrnams.awk: script to correct chromosome names for rmap bed files
# and other files with similar issues
#
# use:
# awk -f <path>tidyrrnams.awk splitfield=<n> <src_file>
#
# splitfield def = 1
#
BEGIN { splitpos=1; splitfield=1; \
#  if (FILENAME=="") \
#    { \
#    printf("use:\n"); \
#    printf("awk splitfield=<n> -f <path>tidyrrnams.awk <src_file>\n"); \
#    printf("  where splitfield is the field to process, def = 1\n"); \
#    } \
  }
{ \
split($splitfield,cnam,"r"); \
  for (i=1; i<= NF; i++) \
    if (i==splitfield) \
      printf("%s%s",((i==1)?"":"\t"),cnam[splitpos]); \
    else \
      printf("%s%s",((i==1)?"":"\t"),$i); \
printf("\n"); \
}


