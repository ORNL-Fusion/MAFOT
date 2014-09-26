// Header-Datei fuer LA_STRING-Klasse

#ifndef LA_STRING_INCLUDED
#define LA_STRING_INCLUDED

//#include "machine.mac"

#include <iostream>
#include <string.h>
#include <ctype.h>

#define LA_STRING_default_case 1    // -> bei Vergleichen Unterscheidung von
							// Gross-/Kleinschreibung
#define LA_STRING_buf_len 100	// Pufferlaenge fuer Eingabe, bei laengeren
							// Eingaben werden Bloecke von jeweils
							// buf_len Zeichen eingelesen
#define LA_STRING_eol '\n'	// Zeilenendezeichen
#define LA_STRING_BOOLEAN int


// Klassendefinition von LA_STRING

class LA_STRING
{
private:

   char *s;                    	// Zeiger auf String (mit abschliessender Null)
   unsigned slen;              	// Laenge des Strings
   static char Buffer[LA_STRING_buf_len];// Puffer fuer Umwandlung/Einlesen

public:

// die folgenden Deklarationen entsprechen der IBM-Klasse istring
// falls nicht anders kommentiert


// Konstruktoren
   inline LA_STRING() {s=new char[1]; s[0]=0; slen=0;};	//  Standardkonstruktor  LA_STRING x
   LA_STRING(const char);         		//  Konstruktor LA_STRING x='a'
   LA_STRING(int);                		//  Konstruktor LA_STRING x=-12
   LA_STRING(unsigned);           		//  Konstruktor LA_STRING x=12
   LA_STRING(const char *const);     		//  Konstruktor LA_STRING x="abc"
   LA_STRING(long);               		//  Konstruktor LA_STRING x=(long)12
   LA_STRING(const LA_STRING&);     		//  Kopierkonstruktor    LA_STRING x=y
   LA_STRING(unsigned long);
   LA_STRING(float, char *f="%g");		//  Konstruktor LA_STRING x=1.2  (nicht IBM)
   LA_STRING(double, char *f="%g");//Konstruktor LA_STRING x=(double)1.2
                                //  (f ist Format-String fuer sprintf)

// Destruktor
   ~LA_STRING() {delete [] s; s=NULL; slen=0;} //  Destruktor


//  Typenumwandlungen
   operator char* () const;          //  typumwandlung in char*

   long   asInt() const;             //  liefert 0, wenn !islong()
   double asDouble() const;          //  liefert 0, wenn !isdouble()

//  Typentests
   LA_STRING_BOOLEAN  isAlphabetic() const; //  liefert 1, wenn String aus A..Z,a..z
   LA_STRING_BOOLEAN  isAlphanumeric() const; //  liefert 1, wenn Alphabetic oder 0..9
//   LA_STRING_BOOLEAN  isASCII() const;      //  liefert 1, wenn ' '..#127
   LA_STRING_BOOLEAN  isDigits() const;     //  liefert 1, wenn 0..9

// Vergleiche

   friend LA_STRING_BOOLEAN operator==(const LA_STRING &s1, const LA_STRING &s2);
   friend LA_STRING_BOOLEAN operator==(const LA_STRING &s1, char s2[]);
   friend LA_STRING_BOOLEAN operator==(char s1[], const LA_STRING &s2);
   friend LA_STRING_BOOLEAN operator==(const LA_STRING &s1, const char *const s2);
   friend LA_STRING_BOOLEAN operator==(const char *const s1, const LA_STRING &s2);

   friend LA_STRING_BOOLEAN operator!=(const LA_STRING &s1, const LA_STRING &s2);
   friend LA_STRING_BOOLEAN operator!=(const LA_STRING &s1, char s2[]);
   friend LA_STRING_BOOLEAN operator!=(char s1[], const LA_STRING &s2);
   friend LA_STRING_BOOLEAN operator!=(const LA_STRING &s1, const char *const s2);
   friend LA_STRING_BOOLEAN operator!=(const char *const s1, const LA_STRING &s2);

   friend LA_STRING_BOOLEAN operator<=(const LA_STRING &s1, const LA_STRING &s2);
   friend LA_STRING_BOOLEAN operator<=(const LA_STRING &s1, const char *const s2);
   friend LA_STRING_BOOLEAN operator<=(const char *const s1, const LA_STRING &s2);
   friend LA_STRING_BOOLEAN operator>=(const LA_STRING &s1, const LA_STRING &s2);
   friend LA_STRING_BOOLEAN operator>=(const LA_STRING &s1, const char *const s2);
   friend LA_STRING_BOOLEAN operator>=(const char *const s1, const LA_STRING &s2);
   friend LA_STRING_BOOLEAN operator<(const LA_STRING &s1, const LA_STRING &s2);
   friend LA_STRING_BOOLEAN operator<(const LA_STRING &s1, const char *const s2);
   friend LA_STRING_BOOLEAN operator<(const char *const s1, const LA_STRING &s2);
   friend LA_STRING_BOOLEAN operator>(const LA_STRING &s1, const LA_STRING &s2);
   friend LA_STRING_BOOLEAN operator>(const LA_STRING &s1, const char *const s2);
   friend LA_STRING_BOOLEAN operator>(const char *const s1, const LA_STRING &s2);
// Achtung Zuweisungen der Form <LA_STRING a="Hi"+'!'> koennen nicht// vorgenommen werden, da die rechte Seite als Zeiger-Addition// uebersetzt wird !

// Operatoren
   friend std::ostream& operator<< (std::ostream&, const LA_STRING&);
          // Ueberladen des Ausgabeoperators
   friend std::istream& operator>> (std::istream&, LA_STRING&);
          // Ueberladen des Eingabeoperators

   LA_STRING &operator =(const LA_STRING&); //  Zuweisungsoperator   x = y
   LA_STRING &operator+=(const LA_STRING&); //  haengt LA_STRING an
//   LA_STRING &operator+=(const char*);
   LA_STRING &operator+=(const char);           //  haengt Zeichen an (nicht IBM)
   friend const LA_STRING operator *(unsigned n,const LA_STRING&);// liefert n-fache Kopie B+..+B (nicht IBM)

   friend const LA_STRING operator +(const LA_STRING&, const LA_STRING&);
   friend const LA_STRING operator +(const LA_STRING&, char[]);
   friend const LA_STRING operator +(char[], const LA_STRING&);
   friend const LA_STRING operator +(const LA_STRING&, const char *const);
   friend const LA_STRING operator +(const char *const, const LA_STRING&);

   friend const LA_STRING operator +(const LA_STRING&, const char);
   friend const LA_STRING operator +(const char, const LA_STRING&);

   const char &operator[] (int i) const; //  Elementzugriff
   char &operator[](int i);              //  Elementzuweisung

   const char &operator[] (unsigned i) const; //  Elementzugriff
   char &operator[](unsigned i);              //  Elementzuweisung

   const LA_STRING operator()(unsigned a,unsigned b) const;
           // liefert b-Zeichen-SuLA_STRING ab a (nicht IBM)


// String-Funktionen
   inline unsigned length() const {return(slen);}; //  Laenge des Strings ohne abschliessende Null

   LA_STRING suLA_STRING( unsigned startPos) const; // Teilstring ab startPos
   LA_STRING suLA_STRING( unsigned startPos, unsigned length, char padChar = ' ') const;
           // Teilstring der Laenge length ab startPos, aufgefuellt mit padChar

   unsigned indexOf(const char, unsigned Start=1) const;
            //  erstes Vorkommen des Zeichens ab Position Start
   unsigned indexOf(const LA_STRING&, unsigned Start=1) const;
            //  erstes Vorkommen des Teilstrings ab Position Start
   unsigned indexOfWord(unsigned n, char Trenner=' ') const;
            // Position des n-ten Wortes  (nicht IBM)

   LA_STRING lineFrom( std::istream& is, char c=LA_STRING_eol); // Zeile aus is einlesen


   LA_STRING_BOOLEAN includes( const LA_STRING &bs) const;
           // Test ob bs als Teilstring vorhanden ist
   LA_STRING_BOOLEAN includes( char ch) const; // Test, ob Zeichen enthalten

   LA_STRING &leftJustify ( unsigned n, char padChar = ' ') const;
           //  liefert die ersten n Zeichen gegebenenfalls Auffuellen
   LA_STRING &rightJustify( unsigned n, char padChar = ' ') const;
           //  liefert die letzten n Zeichen gegebenenfalls Auffuellen
   LA_STRING &strip(char stripChar = ' ');
            // Entfernt alle stripChar am Anfang und Ende des Strings


//  Nicht IBM kompatibel

   LA_STRING_BOOLEAN  isInt() const;     //  testet auf long
   LA_STRING_BOOLEAN  isDouble() const;  //  liefert 1, wenn String ein double darstellt

//  Ein-/Ausgabefunktionen (auch fuer fstream geeignet)
   LA_STRING &readChar(std::istream&); //  Einlesen eines bel. Zeichens
   LA_STRING &readnChar(std::istream&, unsigned); //  Einlesen von genau n Zeichen
   LA_STRING &readLn(std::istream&, char c=LA_STRING_eol); // Einlesen bis Zeilenendezeichen eol
   LA_STRING &readnLn(std::istream&, unsigned, char c=LA_STRING_eol);  // Einlesen von genau n Zeichen oder bis eol
   LA_STRING &nwReadLn(std::istream&, char c=LA_STRING_eol);  // ignoriert Leerzeilen

   inline unsigned len() const {return(slen);};  // Laenge des Strings
   const LA_STRING left  ( unsigned i) const; // liefert die ersten i Zeichen
   const LA_STRING right  ( unsigned i) const;// liefert die letzten i Zeichen
   const LA_STRING mid  ( unsigned i) const;  // liefert restliche Zeichen ab Position i
   const LA_STRING mid  ( unsigned i, unsigned n) const; // liefert n Zeichen von Position i an
   const LA_STRING subst( unsigned i, unsigned n, const LA_STRING& B) const;
           //  ersetzt die n Zeichen ab Position i durch String B
   const LA_STRING insert(const LA_STRING &A, unsigned n=1);
           // fuegt A an Position n ein (Anhaengen: n=len+1)
   const LA_STRING asUpper() const; //  liefert String in Grossbuchstaben
   const LA_STRING asLower() const; //  liefert String in Kleinbuchstaben
   const LA_STRING lowerCase() const; //  liefert String in Kleinbuchstaben
};
// Ende der Definition von LA_STRING



void case_sensitive();		// String-Vergleiche unterscheiden					// Gross-/Kleinschreibung
void not_case_sensitive(); 	// keine Unterscheidung von					// Gross-/Kleinschreibung



#endif // LA_STRING_INCLUDED

