// LA_STRING.cxx
// Implementation der LA_STRING-Klasse

//--------------
// Include-Files
//--------------

#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
/*#include <values.h>
Diese Header-Datei fehlt ??
daher folgende definitionen:
*/
#define MAXINT 32768
#define MAXFLOAT 8.45e+3

#ifndef FALSE
#define FALSE 0
#endif

#ifndef TRUE
#define TRUE 1
#endif


#include "la_string.hxx"

// namespaces
//-----------
using namespace std;


#define error_index "LA_STRING-error: index out of range"
#define error_empty "LA_STRING-error: addressing uninitialized LA_STRING"
#define error_internal "LA_STRING-error: interner Fehler"

int cmp_case=LA_STRING_default_case;                    // Standard-Einstellung fuer Vergleiche
char LA_STRING::Buffer[LA_STRING_buf_len];                // Puffer fuer Umwandlung/Einlesen

void error(const char* s, int LINE)          // Fehlermeldung und Abbruch
{
    cerr << s << '\n';
    cerr << "LA_STRING.cxx: line=" << LINE << '\n';
    abort();
}

#ifndef CHECK
#define CHECK(pc) if(pc==NULL) error(error_empty, __LINE__);
#endif

#ifndef VALIDATE
#define VALIDATE(tst) if(!(tst)) error(error_internal, __LINE__);
#endif

// -------------------------------------------------------------------------------
// Konstruktoren:
// --------------

LA_STRING::LA_STRING(const LA_STRING& x)              // Kopierkonstrunktor   LA_STRING x = y
{  CHECK(x.s)

   slen = x.len();
   s = new char[ slen+1 ];
   strcpy(s,x.s);                             // String kopieren
}


LA_STRING::LA_STRING(const char *const s0)         // Kopierkontruktor     LA_STRING x = "abc"
{
    if (s0==NULL)
    {   s=new char[1]; s[0]=0;
        slen=0;
    }
    else
    {   slen = strlen(s0);
        s = new char[ slen+1 ];
        strcpy(s,s0);                               // String kopieren
    }
}


LA_STRING::LA_STRING(const char c)                  // Kopierkontruktor     LA_STRING x = 'a'
{   if(c)
    {
       s = new char[2];                           // Platz fuer Zeichen schaffen
       slen=1;
       s[0]=c;  s[1]=0;                           // "Zeichen-String" erzeugen
    }
    else
    {
       s = new char[1];                           // Platz fuer Zeichen schaffen
       slen=0;
       s[0]=0;                                    // "Zeichen-String" erzeugen
    }
}


LA_STRING::LA_STRING(int i)
{
    sprintf(Buffer, "%d", i);
    s = new char[ strlen(Buffer)+1 ];
    strcpy(s,Buffer);                                // String kopieren
    slen=strlen(s);
}


LA_STRING::LA_STRING( long i )
{
    sprintf(Buffer, "%ld", i);
    s = new char[ strlen(Buffer)+1 ];
    strcpy(s,Buffer);                                // String kopieren
    slen=strlen(s);
}


LA_STRING::LA_STRING( unsigned i )
{
    sprintf(Buffer, "%u", i);
    s = new char[ strlen(Buffer)+1 ];
    strcpy(s,Buffer);                                // String kopieren
    slen=strlen(s);
}


LA_STRING::LA_STRING(float x)
{
    sprintf(Buffer, "%g", x);
    s = new char[ strlen(Buffer)+1 ];
    strcpy(s,Buffer);                                // String kopieren
    slen=strlen(s);
}

LA_STRING::LA_STRING(double x)
{
    sprintf(Buffer, "%g", x);
    s = new char[ strlen(Buffer)+1 ];
    strcpy(s,Buffer);                                // String kopieren
    slen=strlen(s);
}


//-----------------------------------------------------
// Zuweisungen
// -----------

LA_STRING &LA_STRING::operator=(const LA_STRING &x)   // Zuweisungsoperator   x = y
{
    CHECK(x.s) CHECK(s)

    if(this!=&x)
    {   delete[] s;                        
        slen = x.slen;
        s = new char[ slen+1 ];
        strcpy(s, x.s);                               // String kopieren
    }
    return(*this);
}


LA_STRING &LA_STRING::operator+=(const LA_STRING &A)   // haengt LA_STRING A an
{
    CHECK(A.s) CHECK(s)

    int l1=slen, l2=A.len();
    if (l2==0) return(*this);

    slen = l1+l2;
    char *ss = new char[ slen+1 ];

    if (l1>0)
    {   strcpy( ss, s);                       // erster String nach ss
        strcat( ss, A.s);                     // Summenstring nach ss
    }
    else
        strcpy( ss, A.s);                     // Summenstring nach ss

    if(this!=&A) delete[] s;                  // wird in der Zuweisung erledigt
    s = ss;
    return(*this);
}


LA_STRING &LA_STRING::operator+=(const char A)             // haengt LA_STRING A an
{  CHECK(s)

   if(A==0) return *this;
   if(slen==0) return (*this=A);

   slen++;
   char *ss = new char[ slen+1 ];

   strcpy( ss, s);                       // erster String nach ss
   ss[slen-1] = A;                                  // Summenstring nach ss
   ss[slen] = 0;                                  // String terminieren

   delete[] s; 
   s = ss;
   return (*this);
}

//---------------------------------------------------------------------------
// Member-Funktionen
// -----------------


unsigned LA_STRING::indexOf( const char c, unsigned Start) const
    //  erste Position des Zeichens c ab Position Start
{  CHECK(s)
   if (slen==0) return 0;

   if (Start<1 || Start>slen) error( error_index, __LINE__);
   char *a=strchr( s+Start-1, c);
   if (a!=NULL) return (unsigned)(a-s+1);
   else return(0);
}


unsigned LA_STRING::indexOf(const LA_STRING& B, unsigned Start) const
    //  erste Position eines Teilstrings in s ab Position Start
{  CHECK(s) CHECK(B.s)
   if (Start<1 || Start>slen) error( error_index, __LINE__);
   char *a = strstr( s+Start-1, B.s);
   if (a!=NULL) return (unsigned)(a-s+1);
   else  return(0);
}


unsigned LA_STRING::indexOfWord(unsigned n, char Trenner) const
    // Position des n-ten Worts
{  CHECK(s)
   if(n==0) return(0);
   char *p=s;
   do
   {  while (*p && *p==Trenner) p++;
      if (--n==0) return (unsigned)(p-s+1);
      while (*p && *p!=Trenner) p++;
   } while (*p);
   return 0;
}


const LA_STRING LA_STRING::left ( unsigned n)  const              //  liefert die ersten n Zeichen
{  return mid(1, n);
}


const LA_STRING LA_STRING::right( unsigned n)  const              //  liefert die letzten n Zeichen
{  if (n>=slen) return *this;
   else return mid(slen-n+1);
}

const LA_STRING LA_STRING::mid ( unsigned i) const
{
   if (i<1) error( error_index, __LINE__ );

   return mid(i, slen-i+1);
}

const LA_STRING LA_STRING::mid( unsigned i, unsigned n) const        //  liefert n Zeichen von Position i an
{
   CHECK(s)

   if (i<1) { cerr << "(" << s << ")  Index:" << i << endl << flush; error( error_index, __LINE__ ); }


   LA_STRING B;
   if (n>0)
   {  char *ss = new char[n+1];
      strncpy( ss, s+i-1, n);  ss[n] = 0;
      B = ss;
      delete[] ss;                               // ss wird nicht mehr gebraucht
   }
   return(B);
}


const LA_STRING LA_STRING::subst( unsigned i, unsigned n, const LA_STRING &B) const
{  CHECK(s) CHECK(B.s)
   return(left(i-1)+B+mid(i+n));
}


const LA_STRING LA_STRING::asUpper() const                    //  liefert String in Grossbuchstaben
{  CHECK(s)
   LA_STRING B(s);
   for (unsigned int i=0; i<len(); i++)
     B.s[i]=(char) toupper(B.s[i]);
   return B;
}

const LA_STRING LA_STRING::asLower() const                    //  liefert String in Kleinbuchstaben
{   CHECK(s) 
    LA_STRING B(s);
    for (unsigned int i=0; i<len(); i++)
        B.s[i]=(char) tolower(B.s[i]);
    return B;
}


const LA_STRING LA_STRING::lowerCase() const                    //  liefert String in Kleinbuchstaben
{
    return (*this).asLower();
}


//-----------------------------------------------------------------
// Klassen-Operatoren
// ------------------

const char &LA_STRING::operator[](int i) const           //  Elementzugriff
{
    if (i<1 || i>int(slen) ){ cerr << i; error( error_index, __LINE__ );}
    return( s[i-1] );
}

char &LA_STRING::operator[](int i)   //  Elementzugriff fuer Zuweisung
{
    if (i<1 || i>int(slen) ){ cerr << i; error( error_index, __LINE__ );}
    return( s[i-1] );
}

const char &LA_STRING::operator[](unsigned i) const           //  Elementzugriff
{
    if (i<1 || i>slen ){ cerr << i; error( error_index, __LINE__ );}
    return( s[i-1] );
}

char &LA_STRING::operator[](unsigned i)   //  Elementzugriff fuer Zuweisung
{
    if (i<1 || i>slen ){ cerr << i; error( error_index, __LINE__ );}
    return( s[i-1] );
}

const LA_STRING LA_STRING::operator()(unsigned i, unsigned n) const
{   return mid( i, n);
}


// Operaror ==
LA_STRING_BOOLEAN operator==(const LA_STRING &s1, const LA_STRING &s2) //  testet auf Gleichheit
{  CHECK(s1.s) CHECK(s2.s)
   if (cmp_case) return strcmp( s1.s, s2.s)==0;
            else return strcmp( s1.asUpper().s, s2.asUpper().s)==0;
}

LA_STRING_BOOLEAN operator==(const LA_STRING &s1, char s2[])     //  testet auf Gleichheit
{  CHECK(s1.s)
   if (s2==NULL) return FALSE;

   if (cmp_case)
      return strcmp( s1.s, s2)==0;
   else
    {
      LA_STRING hs(s2);
      return strcmp( s1.asUpper().s, hs.asUpper().s)==0;
    }
}

LA_STRING_BOOLEAN operator==(char s1[], const LA_STRING &s2)
{  return s2==s1;
}

LA_STRING_BOOLEAN operator==(const LA_STRING &s1, const char *const s2)     //  testet auf Gleichheit
{  CHECK(s1.s)
   if (s2==NULL) return FALSE;

   if (cmp_case)
      return strcmp( s1.s, s2)==0;
   else
    {
      LA_STRING hs(s2);
      return strcmp( s1.asUpper().s, hs.asUpper().s)==0;
    }
}

LA_STRING_BOOLEAN operator==(const char  *const s1, const LA_STRING &s2)
{  return s2==s1;
}


// Operator !=
LA_STRING_BOOLEAN operator!=(const LA_STRING &s1, const LA_STRING &s2)
{  return !(s1==s2);
}

LA_STRING_BOOLEAN operator!=(const LA_STRING &s1, char s2[] )
{  return !(s1==s2);
}

LA_STRING_BOOLEAN operator!=(char s1[], const LA_STRING &s2)
{  return !(s2==s1);
}

LA_STRING_BOOLEAN operator!=(const LA_STRING &s1, const char *const s2)
{  return !(s1==s2);
}

LA_STRING_BOOLEAN operator!=(const char *const s1, const LA_STRING &s2)
{  return !(s2==s1);
}


// Operator <=
LA_STRING_BOOLEAN operator<=(const LA_STRING &s1, const LA_STRING &s2)
{ CHECK(s1.s) CHECK(s2.s)
  if (cmp_case) return strcmp( s1.s, s2.s)<=0;
           else return strcmp( s1.asUpper().s, s2.asUpper().s)<=0;
}

LA_STRING_BOOLEAN operator<=(const LA_STRING &s1, const char *const s2)
{  CHECK(s1.s)
   if (s2==NULL) return FALSE;
   if (cmp_case) return strcmp( s1.s, s2)<=0;
   else
   {
    LA_STRING hs(s2);
    return strcmp( s1.asUpper().s, hs.asUpper().s)<=0;
   }
}


LA_STRING_BOOLEAN operator<=(const char *const s1, const LA_STRING &s2)
{  if (s1==NULL) return TRUE;
   return s2>=s1;
}


LA_STRING_BOOLEAN operator>=(const LA_STRING &s1, const LA_STRING &s2) //  testet auf Gleichheit
{  CHECK(s1.s) CHECK(s2.s)
   if (cmp_case) return strcmp( s1.s, s2.s)>=0;
            else return strcmp( s1.asUpper().s, s2.asUpper().s)>=0;
}

LA_STRING_BOOLEAN operator>=(const LA_STRING &s1, const char *const s2)     //  testet auf Gleichheit
{ CHECK(s1.s)
  if(s2==NULL) return TRUE;

  if (cmp_case)
      return strcmp( s1.s, s2)>=0;
  else
    {
      LA_STRING hs(s2);
      return strcmp( s1.asUpper().s, hs.asUpper().s)>=0;
    }
}


LA_STRING_BOOLEAN operator>=(const char *const s1, const LA_STRING &s2)
{  if (s1==NULL) return FALSE;
   return s2<=s1;
}


LA_STRING_BOOLEAN operator<(const LA_STRING &s1, const LA_STRING &s2) //  testet auf Gleichheit
{  return !(s1>=s2);  
}

LA_STRING_BOOLEAN operator<(const LA_STRING &s1, const char *const s2)     //  testet auf Gleichheit
{  return !(s1>=s2);
}

LA_STRING_BOOLEAN operator<(const char *const s1, const LA_STRING &s2)
{  return !(s1>=s2);
}

LA_STRING_BOOLEAN operator>(const LA_STRING &s1, const LA_STRING &s2) //  testet auf Gleichheit
{  return !(s1<=s2);
}

LA_STRING_BOOLEAN operator>(const LA_STRING &s1, const char *const s2)     //  testet auf Gleichheit
{  return !(s1<=s2);
}

LA_STRING_BOOLEAN operator>(const char *const s1, const LA_STRING &s2)
{  return !(s1<=s2);
}


// ---------------------------------------------------------------------------
//  Typenumwandlungen
// ------------------

LA_STRING::operator char* () const
{
   char* p = new char [slen+1];
   strcpy( p, s);
   return(p);
}

long LA_STRING::asInt() const
{  CHECK(s)
   return atol(s);
}

double LA_STRING::asDouble() const                     //  liefert 0, wenn !isfloat()
{  CHECK(s)
   return atof(s);
}

LA_STRING_BOOLEAN  LA_STRING::isInt() const                        //  liefert 1, wenn String ein long darstellt
{  CHECK(s)
   char *err;
   strtol( s, &err, 10);
   return *err==0;
}

LA_STRING_BOOLEAN  LA_STRING::isDouble() const                      //  liefert 1, wenn String ein double darstellt
{  CHECK(s)
   char *err;
   strtod( s, &err);
   return *err==0;
}

LA_STRING_BOOLEAN LA_STRING::isAlphabetic() const            //  liefert 1, wenn String aus A..Z,a..z besteht
{  CHECK(s)
   char *p=s;
   while (*p && isalpha(*p)) p++;
//   while (*p && isascii(*p) && isalpha(*p)) p++;
   return *p == 0;
}

LA_STRING_BOOLEAN LA_STRING::isAlphanumeric() const          //  liefert 1, wenn Alphabetic oder 0..9
{  CHECK(s)
   char *p=s;
   while (*p && isalnum(*p)) p++;
//   while (*p && isascii(*p) && isalnum(*p)) p++;
   return *p == 0;
}

/*
LA_STRING_BOOLEAN LA_STRING::isASCII() const                 //  liefert 1, wenn ' '..#127
{  CHECK(s)
   char *p=s;
   while (*p && isascii(*p)) p++;
   return *p == 0;
}
*/

LA_STRING_BOOLEAN LA_STRING::isDigits() const                //  liefert 1, wenn 0..9
{  CHECK(s)
   char *p=s;
   while (*p && isdigit(*p)) p++;
//   while (*p && isascii(*p) && isdigit(*p)) p++;
   return *p == 0;
}

// ---------------------------------------------------------------------------
//  Ein-/Ausgabefunktionen
//  ----------------------

LA_STRING &LA_STRING::readChar(istream& is)                //  Einlesen eines bel. Zeichens
{
char c;

if (is.get(c))
        *this = c;
   else
        *this = "";
   return *this;
}

LA_STRING &LA_STRING::readnChar(istream& is, unsigned n)
   //  Einlesen von genau n bel. Zeichen
{  CHECK(s)
   if (n<1) return *this;

   delete[] s;
   s = new char[ n+1];
   is.read(s,n);
   VALIDATE((slen=is.gcount())<=n)
   s[slen]=0;
   return *this;
}

LA_STRING &LA_STRING::readLn(istream& is, char eol) //  Einlesen bis Zeilenende
{
*this="";

int p = is.peek();
if (p==eol)
	{
	is.ignore(1);
	p = is.peek();
	}

while ( p !=eol && p !=EOF && !is.eof())
	{
	is.get( Buffer, LA_STRING_buf_len, eol);
	*this += Buffer;
	p = is.peek();
	}
if (p==eol) is.ignore(1);

return (*this);
}


LA_STRING &LA_STRING::readnLn(istream& is, unsigned n, char eol)
   //  Einlesen von genau n Zeichen oder bis eol
{  CHECK(s)
   if (n<1) return *this;

   delete[] s;
   s = new char[ n+1];
   is.get(s,n+1,eol);
   VALIDATE((slen=is.gcount())<=n)
   if (is.peek()==eol) is.get(eol);
   s[slen]=0;
   return *this;
}


LA_STRING &LA_STRING::nwReadLn(istream& is, char eol)
{
  do
     readLn(is, eol);
  while (slen==0 && !is.fail() && !is.eof());
  return(*this);
}

LA_STRING LA_STRING::lineFrom( istream& aStream, char c)
{
  nwReadLn( aStream, c);
  return( *this);
}

////////////////////////////////////////////////////////////////////////////

LA_STRING LA_STRING::suLA_STRING( unsigned startPos) const
{ CHECK(s)
  return(mid(startPos));
}

LA_STRING LA_STRING::suLA_STRING(  unsigned startPos, unsigned length, char padChar) const
{  CHECK(s)
   LA_STRING bs(mid(startPos,length));
   if (bs.len()<length)
   {
     bs=bs+(length-bs.len())*LA_STRING(padChar);
   }
   return( bs);
}


LA_STRING_BOOLEAN LA_STRING::includes( const LA_STRING  &bs) const
{  CHECK(s)
   LA_STRING_BOOLEAN b = indexOf( bs);
   return (b!=0);
}


LA_STRING& LA_STRING::strip(char stripChar)
{  CHECK(s)
   int L(0),R(slen-1);
   while(L<int(slen) && s[L]==stripChar) L++;                      
   while(R>=0 && s[R]==stripChar) R--;
   if(R<L) *this="";
   for(int i=0; i+L<=R; i++) s[i]=s[i+L];
   s[R-L+1]=0;
   slen = R-L+1;
   return *this;
}


ostream& operator<< (ostream &os, const LA_STRING &B)     // Ueberladen des Ausgabeoperators
{  CHECK(B.s)
   return os << B.s;
}

istream& operator>> (istream &is, LA_STRING &B)           // Ueberladen des Eingabeoperators
{  CHECK(B.s)
   B.readLn(is);
   return(is);
}

const LA_STRING operator* (unsigned n, const LA_STRING& B)           // liefert n-fache Kopie B+..+B
{  CHECK(B.s)
   LA_STRING A;
   for (unsigned int i=0; i<n; i++)  A+=B;
   return A;
}

// operator +
const LA_STRING operator+ (const LA_STRING& s1,const LA_STRING& s2)   // Summenstring
{ CHECK(s1.s) CHECK(s2.s)
  LA_STRING S(s1);
  S+=s2;
return S;
}

const LA_STRING operator+ (const LA_STRING& s1, char p[])
{
   CHECK(s1.s)
   if(p==NULL) 
      return s1;
   else
   {
     LA_STRING bs(s1);
     bs += LA_STRING(p);
     return(bs);
   }
}

const LA_STRING operator+ ( char p[], const LA_STRING& s2 )
{  CHECK(s2.s)
   if(p==NULL) 
      return s2;
   else
   {
     LA_STRING bs(p);
     bs += s2;
     return(bs);
   }
}

const LA_STRING operator+ (const LA_STRING& s1, const char *const p)
{
   CHECK(s1.s)
   if(p==NULL) 
      return s1;
   else
   {
     LA_STRING bs(s1);
     bs += LA_STRING(p);
     return(bs);
   }
}

const LA_STRING operator+ ( const char *const p, const LA_STRING& s2 )
{  CHECK(s2.s)
   if(p==NULL) 
      return s2;
   else
   {
     LA_STRING bs(p);
     bs += s2;
     return(bs);
   }
}

const LA_STRING operator+ (const LA_STRING& bs1, const char c)
{  CHECK(bs1.s)
   LA_STRING bs(bs1);
   bs += c;
   return(bs);
}

const LA_STRING operator+ ( const char c,const LA_STRING& bs2 )
{  CHECK(bs2.s)
   LA_STRING bs(c);
   bs += bs2;
   return(bs);
}

const LA_STRING LA_STRING::insert(const LA_STRING &A, unsigned n) 
//  liefert s mit A an Position n
{  CHECK(A.s) CHECK(s)
   return *this=subst(n,0,A);
}


void case_sensitive()                // String-Vergleiche unterscheiden Gross-/Kleinschreibung
{  cmp_case = 1;  }

void not_case_sensitive()            // keine Unterscheidung von Gross-/Kleinschreibung
{  cmp_case = 0;  }


