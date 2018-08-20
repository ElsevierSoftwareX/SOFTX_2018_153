
Copyright (c) 2018 Gabriel Johnson <gsjohnso@go.olemiss.edu>. All rights
reserved.

1. Introduction
2. Setting up
3. Running
4. Using custom paths to Mathematica and Clozure
5. Using the source code as a library


1. Introduction

This is a program written for mathematical research--specifically that of Dr.
Sandra Spiroff at the University of Mississippi--that calculates the
Hilbert-Kunz multiplicity and F-signature for a given intersection algebra.

This is something of a maiden voyage for me, so I apologize for failing to
follow usual coding practices in some respects. In particular, there's a chance
you might have to edit some of the scripts in it get it to work. I'll be working
on polishing it more in future versions. Sorry for the inconvenience.


2. Setting up

This program assumes that you have either Mac OS X, Linux, or Windows with
Cygwin or Git Bash, as well as Mathematica. If you don't have Mathematica, then
you can still use the program to generate the Mathematica code that would
calculate the Hilbert-Kunz multiplicity and F-signature, but you won't be able
to actually calculate concrete values.

To start, you need to have Clozure Common Lisp (the interpreter for the
programming language used in this project) in the program root directory.

Unless you downloaded the source code directly from Gitlab, it should already
come with Clozure, and you can skip the next few paragraphs. If you did get it
from Gitlab, then read on to see how to install Clozure, or tell the program
where to locate an existing installation.

If you are a programmer with very good taste in programming languages and
already have Clozure installed somewhere else, then you can replace the
ccl-name script under the src directory with your own script that outputs the
path to your installation, as described below in section 4.

If not, you can get Clozure Common Lisp from:
https://ccl.clozure.com/download.html

When you extract the files from the tarball/zip file, the result might be either
a directory called something like "ccl-1.11.5-linuxx86" (possibly with a
different version number and/or OS name), or just a directory called "ccl". What
you want is for the directory called just "ccl" to be in the root directory of
this program. If you got the longer version, then "ccl" should be inside that
directory.

Once you have Clozure in place, you can compile the code from the terminal by
running the script "setup" in the program root directory. To do this, open a
terminal, cd into the program's root directory, and enter

./setup

After running the script, a file called "hkm.image" should appear. You'll need
this file to run the main programs, calculate-integral and inequalities. If
hkm.image were ever to get lost for some reason, then you'll have to run "setup"
again. The setup script has no side effects other than generating hkm.image, so
there's no need to be afraid of running it multiple times.

You can't run the main program without running "setup" first. If you try to, you
should get an error saying something to the effect of, "hkm.image could not be
found."


3. Running

To run the program, open a terminal, cd to the program root directory, and enter

./calculate-integral [a1 a2 ... an] [b1 b2 ... bn]

where each a1 ... an and b1 ... bn should be a positive integer. The numbers
should be separated by spaces. The numbers will be divided into two lists, but
there is no need to give any kind of indicator of the boundary between them; the
list of all the inputs will be split down the middle automatically. Each of
these lists is taken respectively to be the exponents of the generators of an
ideal of a polynomial ring.

The exponents should be in fan order. This means that the quotients a1/b1,
a2/b2, ..., an/bn should be ordered from least to greatest. Geometrically
speaking, if you take each pair of exponents as a point in the first quadrant of
the coordinate plane, as (b1, a1), (b2, a2), ..., (bn, an), then they should all
be ordered counterclockwise starting from the positive x axis.

If any of these conditions is violated, then the program should terminate with
an error message. Otherwise, it outputs the Hilbert-Kunz multiplicity and
F-signature of the intersection algebra of the two ideals.

The program uses Mathematica to calculate the final values of the Hilbert-Kunz
multiplicity and F-signature. However, we've found that starting at around six
or eight inputs, the amount of time Mathematica takes to calculate the values
varies wildly depending on how naturally complex the given intersection algebra
is, often taking as long as several hours on a normal computer. If you are not
willing to wait that long (or risk your computer freezing), or if you just don't
have Mathematica, then you can see the Mathematica code without calculating the
final value by entering

./inequalities [a1 a2 ... an] [b1 b2 ... bn]

The same rules as above apply for what inputs are valid.

With this, you can, for example, run the code on a supercomputer to obtain
results that you wouldn't be able to on a normal computer, or 3D print the
solids whose volumes give the Hilbert-Kunz multiplicity and F-signature.

The inequality-generating part of the program has finished quickly for all
inputs tested so far, so this should pretty much always be safe.

Examples:

# --- Correct usages ---

$ ./calculate-integral 3 2   # the intersection algebra of (x^3) and (x^2)
Hilbert-Kunz Multiplicity = 41/18
F-Signature               = 11/36

$ ./calculate-integral 5 2 2 3   # the intersection algebra of (x^5*y^2)
                                 # and (x^2*y^3)
Hilbert-Kunz Multiplicity = 37283/9900
F-Signature               = 1087/29700

$ ./inequalities 3 2   # without the final integral; the first expression is for
                       # the Hilbert-Kunz multiplicity, the second is for the
                       # F-signature.
Integrate[Boole[
0 <= x && 0 <= y && 3x <= z1 && 2y <= z1 && 
((z1 < 3x + 1 || z1 < 2y + 1)) &&
((x < 1 || y < 2 || z1 < 3x + 1)) &&
((y < 1 || z1 < 3x + 2)) &&
((x < 1 || z1 < 2y + 3)) &&
((x < 1 || y < 1 || z1 < 2y + 1)) &&
((x < 2 || y < 3))], {x, 0, 2000}, {y, 0, 2000}, {z1, 0, 2000}]

Integrate[Boole[
3x <= z1 <= 1 + 3x &&
2y <= z1 <= 1 + 2y], {x, 0, 1}, {y, 0, 1}, 
{z1, 0, 30}]

# --- Incorrect usages ---

$ ./calculate-integral 2 3 5 2  # not in fan order because 5/2 > 2/3
ccl/lx86cl64: Failed with error: Provided numbers not in fan order.

$ ./calculate-integral 5 2 3
ccl/lx86cl64: Failed with error: Odd number of arguments.


4. Using custom paths to Mathematica or Clozure

The scripts that run the program assume that Mathematica is in a standard
location depending on your OS. If you installed Mathematica in a custom
location, then the script won't find it. Likewise, the scripts also assume that
Clozure is in the root directory of the program and will not detect it if you
already have it installed somewhere else.

Probably the easiest way to get around this is to replace the scripts
wolfram-path and ccl-name.

For Mathematica, just go to the src directory, delete wolfram-path, and create a
new text file containing the text below between the lines of hyphens (not
including the hyphens themselves):

--------------------------------------------------------------------------------
#!/bin/bash

echo "/path/to/wolfram"
--------------------------------------------------------------------------------

Replace the text between the quotes with the path to the executable file called
"wolfram" on your machine ("MathKernel" should also work, unless you're on
Windows). Then save it to the src directory under the name "wolfram-path", and
enter the following command at a terminal in the src directory:

chmod +x wolfram-path

To see if it worked, enter

./wolfram-path

and make sure it outputs the path to Mathematica.

Likewise, to use your own path to Clozure, go to the src directory, delete
ccl-name, and replace it with a new text file with the name "ccl-name". Fill in
the contents according to the same template as above, but replace the text in
quotes with the path to Clozure instead of Mathematica. Then enter

chmod +x ccl-name

at a terminal in the src directory. To make sure it worked, enter

./ccl-name

and make sure it outputs the path to Clozure.

NOTE FOR WINDOWS USERS: Whereas text files on Mac OS X and Linux use a special
character called a line feed to indicate the end of a line, Windows uses a
combination of two characters, a carriage return and a line feed. The trouble
is, if you try to run a script that uses Windows-style line endings on Cygwin,
Cygwin will see the extra carriage return character, not realize that it's meant
to indicate a new line, and interpret it as a "word" of its own. If you replace
the scripts in this project, try to run them, and get an error message like,

Command not found: $'\r'

then this is Cygwin tripping over the carriage returns. To get the scripts to
work, you'll have to convert them to Linux-/Mac-style line endings.

One utility for doing this is a program called dos2unix. Alternatively, you
could install a Linux text editor like Nano, Vim, or Emacs into Cygwin (if you
haven't already) by re-running the Cygwin installer. Of these three, Nano is
probably the most beginner-friendly.


5. Using the code as a library

The Lisp portion of the program is complete on its own and can be used as a
library.

But this is where I made what was probably my biggest breach of standard
practice: I didn't put the code in a separate package, create an ASDF system, or
use conditional compilation or libraries to make the code
implementation-independent.

My reasons for doing this were:

- The primary audience for this program is extremely niche, to say the least; it
is not likely to attract a large number of programmers looking to work on it.
Further, it is likely that many of the users of this program are not very
computer savvy. Hence I priorized being able to give the user very specific
instructions on how to set up and run the program rather than giving many
options of implementations.

- Using ASDF and Quicklisp to download external libraries seemed like more of a
cost than a benefit. Considering the audience of this project, the user likely
would not have their own fully furnished Lisp installation to begin with. Thus
on top of Lisp itself, I would need to have them install ASDF and Quicklisp as
well. This just adds extra steps and takes up extra disk space for a user who
would likely never use them for anything else. It also makes the coding more
complicated for me, as I would have to learn a new library's methods of
retrieving command line arguments, exiting the program, and generating image
files rather than using the functions that come with my implementation. Also,
for the compilation process, I would have to write code not just to generate the
image file, but also to make sure that ASDF and Quicklisp are installed.
Foregoing ASDF and Quicklisp spares the user an extra download and makes the
coding simpler.

- Considering that all the code fits into a single file, creating an ASDF system
seems like overkill.

So for these reasons, I took the easier route of making the code
implementation-dependent and not using ASDF or Quicklisp. I apologize to anyone
who might find this inconvenient.

Probably the easiest thing to do right now would be to copy the source code file
into your own project (possibly deleting the function #'main, since it uses
symbols from the CCL package) and just adding a defpackage and in-package form
to the top of the file.
