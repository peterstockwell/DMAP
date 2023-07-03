/* dblu_fns.c: seq database lookup functions in c */

#include <stdlib.h>
#include <stdio.h>
#include <strings.h>
#include <string.h>
/* #include <sys/mode.h> */
#include <sys/types.h>
#include <sys/stat.h>
#include <gdbm.h>
#include <fcntl.h>

#include "bas_fns.h"
#include "wlu_fns.h"
#include "sqfl_fns.h"
#include "sqtrans.h"
#include "sqmat_fns.h"
#include "dbpars.h"
#include "dblu_fns.h"

/* external data for tracking gdbm errors */

extern gdbm_error gdbm_errno;

DB_FLSTREC *db_appnfrec(DB_FLSTREC **spt,
                        char *nm,
                        int no,
                        DBLU_DBFMT fmt)
/* malloc and append a new record to list *spt.  return the address of it */
{
DB_FLSTREC *prev, *endp;

if (spt != NULL)
  {                     /* chain to end of list */
  prev = endp = *spt;
  while (endp != NULL)
    {
    prev = endp;
    endp = endp->nxtfrec;
    }
  endp = (DB_FLSTREC *) getmemory(sizeof(DB_FLSTREC),"File record element");
  endp->nxtfrec = NULL;
  endp->fnam = nm;
  endp->fno = no;
  endp->ffmt = fmt;
  if (*spt == NULL)
    {
    *spt = endp;
    endp->prvfrec = NULL;
    }
  else
    {
    prev->nxtfrec = endp;
    endp->prvfrec = prev;
    }
  }
return(endp);
}

void db_delfrec(DB_FLSTREC *sp,
                DB_FLSTREC **lstrt)
/* delete sp from list *lstrt */
{
DB_FLSTREC *pt;

if (sp != NULL)
  {
  if ((pt = sp->prvfrec) == NULL)
    *lstrt = sp->nxtfrec;
  else
    pt->nxtfrec = sp->prvfrec;
  if ((pt = sp->nxtfrec) != NULL)
    pt->prvfrec = sp->prvfrec;
  memfree(sp->fnam);
  memfree(sp);
  }
}

void db_killflst(DB_FLSTREC **flst)
  /* iteratively delete *flst */
{
while (*flst != NULL)
  db_delfrec(*flst,flst);
*flst = NULL;
}

DB_FLSTREC *db_flstrec4int(DB_FLSTREC *flst,
                           int no)
/* return pointer to the record matching no, NULL if none */
{
DB_FLSTREC *fp;

fp = flst;
while (fp != NULL)
  if (fp->fno == no)
    return(fp);
  else
    fp = fp->nxtfrec;
return(NULL);
}

int db_maxno4flst(DB_FLSTREC *flst)
  /* interate thru flst, return largest no present */
{
int mx;
DB_FLSTREC *fp;

fp = flst;
mx = 0;
while (fp != NULL)
  {
  if (fp->fno > mx)
    mx = fp->fno;
  fp = fp->nxtfrec;
  }
return(mx);
}

size_t db_flkeyhdr(unsigned int fkey,
                   char *ustr,
                   int ulen)
/* generate the string for a file ndbm entry, for fkey. write to ustr up to
limit ulen.  Return no of bytes in total */
{
char *sp;
char *dp;
int acnt;

dp = ustr;
sp = "__DB**NAM";
while (*sp != '\0')
  bas_appchr(ustr,&dp,*(sp++),ulen);
acnt = 0;
while ((fkey > 0) || (acnt == 0))
  {
  bas_appchr(ustr,&dp,(char) (fkey % 0xff) ,ulen);
  fkey /= 0xff;
  acnt++;
  }
return(dp - ustr);
}

size_t db_dbvershdr(char *ustr,
                    int ulen)
/* generate the string for a dbm version entry, for fkey. write to ustr up to
limit ulen.  Return no of bytes in total */
{
char *sp;
char *dp;

dp = ustr;
sp = "__DB**VERS";
while (*sp != '\0')
  bas_appchr(ustr,&dp,*(sp++),ulen);
return(dp - ustr);
}

DB_FLSTREC *db_chkindxdfls(GDBM_FILE dbh)
  /* look specifically for file records and store in linked list: return
the address of the list */
{
datum kdat;
unsigned int fcnt;
datum ddat;
char *kstr;
DB_FLSTREC *flst;
DBLU_DBFMT fmt;
DB_FLSTREC *fend;

fcnt = 1;
fend = flst = NULL;
kstr = (char *) getmemory((DB_MAXKEYLN+DB_FNODIGITS+1),"ndbm key string");
do
  {
  kdat.dptr = kstr;
  kdat.dsize = db_flkeyhdr(fcnt,kdat.dptr,DB_MAXKEYLN+DB_FNODIGITS);
  ddat = gdbm_fetch(dbh,kdat);
  if (ddat.dptr == NULL)
    fcnt = 0;
  else
    {
    fmt = db_chr2dbfmt(*(ddat.dptr));
    fend = db_appnfrec(&flst,db_strdupn(ddat.dptr+2,ddat.dsize-2),fcnt,fmt);
    fcnt++;
    }    
  }
while (fcnt > 0);
free(kstr);
return(flst);
}

void db_saygdbmerr(const char *str)
  /* to pass to gdbm_open as fatal_func() */
{
fprintf(stderr,"gdbm err: %s\n",((str != NULL)?str:"NULL"));
}

DB_NDXTYPE *db_opndbmdb(char *nm,
                        int ro,
                        int oldstyl)
/* open and return the data structure handle for a gdbm database with the 
name handle nm.  if ro (=readonly), then open in an appropriate mode.
oldstyl forces the older form of index   */
{
GDBM_FILE db;
DB_NDXTYPE *dbs;
datum vdat;
datum hdat;
char *kstr;
int prexsts;
char *vstr;

prexsts = (db = gdbm_open(nm,512,GDBM_READER,0,db_saygdbmerr)) != NULL;
if (!ro)   /* must close and reopen */
  {
  if (db != NULL)
    gdbm_close(db);
  db = gdbm_open(nm,512,GDBM_NEWDB,(S_IRUSR|S_IWUSR|S_IRGRP|S_IROTH),
                    db_saygdbmerr);
  }
if (db != NULL)
  {
  dbs = (DB_NDXTYPE *) getmemory(sizeof(DB_NDXTYPE),"gdb file structure");
  dbs->dbh = db;
  dbs->fllst = NULL;
  dbs->locbyt = 4;
  dbs->flbyt = 2;
  hdat.dptr = kstr = (char *) getmemory((DB_MAXKEYLN+1),
                                           "ndbm version string");
  hdat.dsize = db_dbvershdr(hdat.dptr,DB_MAXKEYLN);
  if (prexsts)
    {
    vdat = gdbm_fetch(db,hdat);
    if (vdat.dptr == NULL)        /* no version key, old file */
      {
      dbs->fvers = 0.0;
      dbs->locbyt = 4;
      dbs->flbyt = 1;
      }
    else
      dbs->fvers = (float) atof(vdat.dptr);
    dbs->fllst = db_chkindxdfls(db);
    }
  else
    if (oldstyl)        /* just continue with old format */
      {
      dbs->fvers = 0.0;
      dbs->locbyt = 4;
      dbs->flbyt = 1;
      }
    else
      {      /* new db: assert version no in it */
      dbs->fvers = DBLUVERSIONNO;
      vdat.dptr = vstr = (char *) getmemory((DB_MAXKEYLN+1),
                                             "ndbm version string");
      sprintf(vdat.dptr,"%.2f",DBLUVERSIONNO);
      vdat.dsize = strlen(vdat.dptr);
      if ((gdbm_store(db,hdat,vdat,GDBM_REPLACE)) != 0)
        fprintf(stderr,"Version info store failed for %s\n",nm);
      memfree(vstr);
      }
  memfree(kstr);
  return(dbs);
  }
else
  {
  fprintf(stderr,"GDBM open failed for %s: reason=%d(%s)\n",nm,
            (int) gdbm_errno,gdbm_strerror(gdbm_errno));
  return(NULL);
  }
}

void db_closegdbmdb(DB_NDXTYPE *dbt)
  /* close gdbm database, relinquish storage for dbt */
{
gdbm_close(dbt->dbh);
db_killflst(&dbt->fllst);
dbt->dbh = NULL;
dbt->fllst = NULL;
}

unsigned long db_str2offst(char *str,
                           int sln)
/* scan str for offset value */
{
unsigned long rv;
char *sp;

rv = 0;
sp = str;
while (sln > 0)
  {
  rv = rv << 8;
  rv += (unsigned char) *sp;
  sp++;
  sln--;
  }
return(rv);
}
                  
DBT_DBMTYPE dbt_classdbent(char *dstr)
/* identify dstr as file/id or accno */
{
char *kstr;
int ksln;
DBT_DBMTYPE rvl;

rvl = DBT_db_id;
kstr = (char *) getmemory((DB_MAXKEYLN+DB_FNODIGITS+1),"gdbm file key string");
ksln = db_flkeyhdr(0,kstr,DB_MAXKEYLN+DB_FNODIGITS) - 1;
if (bas_strncmp(kstr,dstr,ksln) == 0)
  rvl = DBT_db_file;
else
  {
  ksln = db_dbvershdr(kstr,DB_MAXKEYLN);
  if (bas_strncmp(kstr,dstr,ksln) == 0)
    rvl = DBT_db_versn;
  else
    if (*dstr == '_')
      rvl = DBT_db_acc;
  }
memfree(kstr);
return(rvl);
}

int dbt_showdbents(FILE *ofl,
                   DB_NDXTYPE *dbnt,
                   int lln,
                   DBL_LSTMODE lstyle)
/* show all entries non-file entries from dbh at ofl, according to lstyle.
Return number shown */
{
datum kdat;
int fcnt;
datum ddat;
DBLU_DBFMT fmt;
DB_FLSTREC *fend;
DBT_DBMTYPE dbt;
int cc;
int acnt;

cc = fcnt = 0;
kdat = gdbm_firstkey(dbnt->dbh);
switch (lstyle)
  {
  case DBL_lst_full:
    while (kdat.dptr != NULL)
      {
      ddat = gdbm_fetch(dbnt->dbh,kdat);
      if (ddat.dptr != NULL)
        switch (dbt = dbt_classdbent(kdat.dptr))
          {
          case DBT_db_id:
          case DBT_db_acc:
            if (dbt == DBT_db_id)
              fputs("ID: ",ofl);
            else
              fputs("ACC: ",ofl);
            fcnt++;
            (void) bas_coutputstr(ofl,kdat.dptr,'"',kdat.dsize,bas_coutputchr);
            fputc(' ',ofl);
            (void) bas_coutputstr(ofl,ddat.dptr,'"',ddat.dsize,bas_putochr);
            fprintf(ofl,"(=File %d at %lu)\n",
                      (int) db_str2offst(ddat.dptr,dbnt->flbyt),
                      db_str2offst((ddat.dptr+dbnt->flbyt),
                                     (ddat.dsize-dbnt->flbyt)));
            break;
          case DBT_db_versn:
            fprintf(ofl,"VERS: %s\n",ddat.dptr);
            break;
          case DBT_db_file:
          default:
            break;
          }
      kdat = gdbm_nextkey(dbnt->dbh,kdat);
      }
    break;
  case DBL_lst_idlst:
    while (kdat.dptr != NULL)
      {
      if (dbt_classdbent(kdat.dptr) == DBT_db_id)
        {
        ddat = gdbm_fetch(dbnt->dbh,kdat);
        if (ddat.dptr != NULL)
          {
          dbt_putidstr2fl(ofl,lln,&cc,kdat.dptr,kdat.dsize);
          fputc('\n',ofl);
          cc = 0;
          fcnt++;
          }
        }
      kdat = gdbm_nextkey(dbnt->dbh,kdat);
      }
    break;
  case DBL_lst_norm:
  default:
    while (kdat.dptr != NULL)
      {
      if (dbt_classdbent(kdat.dptr) == DBT_db_id)
        {
        ddat = gdbm_fetch(dbnt->dbh,kdat);
        if (ddat.dptr != NULL)
          {
          if (fcnt <= 0)
            dbt_putstr2fl(ofl,lln,&cc,"Id:");
          dbt_putidstr2fl(ofl,lln,&cc,kdat.dptr,kdat.dsize);
          fcnt++;
          }
        }
      kdat = gdbm_nextkey(dbnt->dbh,kdat);
      }
    if (cc > 0)
      fputc('\n',ofl);
    acnt = cc = 0;
    kdat = gdbm_firstkey(dbnt->dbh);
    while (kdat.dptr != NULL)
      {
      if (dbt_classdbent(kdat.dptr) == DBT_db_acc)
        {
        ddat = gdbm_fetch(dbnt->dbh,kdat);
        if (ddat.dptr != NULL)
          {
          if (acnt <= 0)
            dbt_putstr2fl(ofl,lln,&cc,"Acc:");
          dbt_putidstr2fl(ofl,lln,&cc,(kdat.dptr+dbnt->flbyt),
                                         (kdat.dsize-dbnt->flbyt));
          acnt++;
          fcnt++;
          }
        }
      kdat = gdbm_nextkey(dbnt->dbh,kdat);
      }
    if (cc > 0)
      fputc('\n',ofl);
    break;
  }
return(fcnt);
}

DB_FLSTREC *db_dbminfo4nm(DB_NDXTYPE *dbn,
                          char *ustr,
                          datum *edat)
/* return pointer to flst element for ustr corresponding to id/accno ustr.
edat is used for the output datum.  Return NULL for not found. */
{
datum kdat;
int fno;

kdat.dptr = ustr;
kdat.dsize = strlen(ustr);
*edat = gdbm_fetch(dbn->dbh,kdat);
if (edat->dptr == NULL)
  return(NULL);
else
  {
  fno = (int) db_str2offst(edat->dptr,dbn->flbyt);
  return(db_flstrec4int(dbn->fllst,fno));
  }
}

int db_offst2str(unsigned long offst,
                 char *str,
                 int sln)
/* attempt to encode offst into str in sln chars.  return 1 if succeeds */
{
char *sp;

sp = str + sln - 1;
while (sp >= str)
  {
  *sp = (unsigned char) 0xff & offst;
  offst = offst >> 8;
  sp--;
  }
return(offst == 0);
}

FILE *db_opnnsetfile4frp(DB_NDXTYPE *dbn,
                         DB_FLSTREC *frp,
                         datum *edat,
                         DBLU_DBFMT *fmt,
                         char **fnam)
/* given frp & edat, attempt to open the file and position to correct point.
Return NULL for failure.
Set fmt & *fnam from frp */
{
long offst;
FILE *fret;

if (frp == NULL)
  return(NULL);
else
  if ((fret = fopen(frp->fnam,"r")) != NULL)
    {
    *fmt = frp->ffmt;
    *fnam = frp->fnam;
    if (fseek(fret,db_str2offst((edat->dptr+dbn->flbyt),
                                    (edat->dsize-dbn->flbyt)),SEEK_SET) == 0)
      return(fret);
    else
      {
      fclose(fret);
      return(NULL);
      }
    }
  else
    return(NULL);
}

FILE *db_setfile4unam(DB_NDXTYPE *dbn,
                      char *ustr,
                      DBLU_DBFMT *fmt,
                      char **fnam)
/* if ustr can be located as a datum in dbn->dbh, then attempt to open the
corresponding file at the appropriate position. Establish the format of it 
and return the file handle.
Return NULL for failure */
{
datum ddat;
DB_FLSTREC *frp;

if ((dbn == NULL) || (dbn->dbh == NULL))
  return(NULL);
else
  if ((frp = db_dbminfo4nm(dbn,ustr,&ddat)) == NULL)
    return(NULL);
  else
    return(db_opnnsetfile4frp(dbn,frp,&ddat,fmt,fnam));
}

void dbt_showflhdrs(FILE *ofl,
                    DB_FLSTREC *flst)
/* report file lines in flst to ofl */
{
DB_FLSTREC *fp;

fp = flst;
while (fp != NULL)
  {
  fprintf(ofl,"%d %s ",fp->fno,db_dbfmt2str(fp->ffmt));
  bas_coutputstr(ofl,fp->fnam,'"',strlen(fp->fnam),bas_coutputchr);
  fputc('\n',ofl);
  fp = fp->nxtfrec;
  }
}
