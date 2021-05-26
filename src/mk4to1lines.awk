# mk4to1lines.awk: take standard in and write out each 4 lines concated
# into 1, with catchar (defaults to "~") separating each concatted line
BEGIN { if (catchar == "") catchar = "~"; lcnt=1 }
{if (lcnt <= 3)
  {
  printf("%s%s",$0,catchar);
  lcnt++;
  }
else
  {
  printf("%s\n",$0);
  lcnt = 1;
  }
}
END { if (lcnt > 1) printf("\n"); }
