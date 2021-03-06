# celltocellreproduce
Code by Brian Camley (https://bcamley.github.io)
README - 

This package contains MATLAB code and data that should be sufficient to reproduce everything in the paper, "Cell-to-cell variation sets a tissue-rheology-dependent bound on collective gradient sensing", currently available on the arxiv at: https://arxiv.org/abs/1707.03532


Each folder provides that code and data for a different element of the paper:

reproduce_mle_clean: the Maximum Likelihood Estimation results from Fig. 1\
reproduce_rot_clean: finding the optimal rotation speed, etc. (Fig. 2)\
reproduce_del_clean: Self-propelled particle simulations using the Delaunay triangulation interactions (Fig. 4, Fig.5)

Within the first two folders, there is a script named reproduce_(something).m. Run this script, and it will prompt you for a variable rerun; if rerun = 1, this script will recompute everything. If rerun = 0, the script will just replot the figures using the pre-generated data.

The third part is more complex and has its own README.

This code is released under the spirit of Matt Might's CRAPL (http://matt.might.net/articles/crapl/) - that code is generally better out than in, even though it may be quite ugly.

                        THE CRAPL v0 BETA 1


0. Information about the CRAPL

If you have questions or concerns about the CRAPL, or you need more
information about this license, please contact:

   Matthew Might
   http://matt.might.net/


I. Preamble

Science thrives on openness.

In modern science, it is often infeasible to replicate claims without
access to the software underlying those claims.

Let's all be honest: when scientists write code, aesthetics and
software engineering principles take a back seat to having running,
working code before a deadline.

So, let's release the ugly.  And, let's be proud of that.


II. Definitions

1. "This License" refers to version 0 beta 1 of the Community
    Research and Academic Programming License (the CRAPL). 

2. "The Program" refers to the medley of source code, shell scripts,
    executables, objects, libraries and build files supplied to You,
    or these files as modified by You.

   [Any appearance of design in the Program is purely coincidental and
    should not in any way be mistaken for evidence of thoughtful
    software construction.]

3. "You" refers to the person or persons brave and daft enough to use
    the Program.

4. "The Documentation" refers to the Program.

5. "The Author" probably refers to the caffeine-addled graduate
    student that got the Program to work moments before a submission
    deadline.


III. Terms

1. By reading this sentence, You have agreed to the terms and
   conditions of this License.
  
2. If the Program shows any evidence of having been properly tested
   or verified, You will disregard this evidence.

3. You agree to hold the Author free from shame, embarrassment or
   ridicule for any hacks, kludges or leaps of faith found within the
   Program.

4. You recognize that any request for support for the Program will be
   discarded with extreme prejudice.

5. The Author reserves all rights to the Program, except for any
   rights granted under any additional licenses attached to the
   Program.


IV. Permissions

1. You are permitted to use the Program to validate published
   scientific claims.

2. You are permitted to use the Program to validate scientific claims
   submitted for peer review, under the condition that You keep
   modifications to the Program confidential until those claims have
   been published.
 
3. You are permitted to use and/or modify the Program for the
   validation of novel scientific claims if You make a good-faith
   attempt to notify the Author of Your work and Your claims prior to
   submission for publication.
 
4. If You publicly release any claims or data that were supported or
   generated by the Program or a modification thereof, in whole or in
   part, You will release any inputs supplied to the Program and any
   modifications You made to the Progam.  This License will be in
   effect for the modified program.


V. Disclaimer of Warranty

THERE IS NO WARRANTY FOR THE PROGRAM, TO THE EXTENT PERMITTED BY
APPLICABLE LAW. EXCEPT WHEN OTHERWISE STATED IN WRITING THE COPYRIGHT
HOLDERS AND/OR OTHER PARTIES PROVIDE THE PROGRAM "AS IS" WITHOUT
WARRANTY OF ANY KIND, EITHER EXPRESSED OR IMPLIED, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE. THE ENTIRE RISK AS TO THE QUALITY AND
PERFORMANCE OF THE PROGRAM IS WITH YOU. SHOULD THE PROGRAM PROVE
DEFECTIVE, YOU ASSUME THE COST OF ALL NECESSARY SERVICING, REPAIR OR
CORRECTION.


VI. Limitation of Liability

IN NO EVENT UNLESS REQUIRED BY APPLICABLE LAW OR AGREED TO IN WRITING
WILL ANY COPYRIGHT HOLDER, OR ANY OTHER PARTY WHO MODIFIES AND/OR
CONVEYS THE PROGRAM AS PERMITTED ABOVE, BE LIABLE TO YOU FOR DAMAGES,
INCLUDING ANY GENERAL, SPECIAL, INCIDENTAL OR CONSEQUENTIAL DAMAGES
ARISING OUT OF THE USE OR INABILITY TO USE THE PROGRAM (INCLUDING BUT
NOT LIMITED TO LOSS OF DATA OR DATA BEING RENDERED INACCURATE OR
LOSSES SUSTAINED BY YOU OR THIRD PARTIES OR A FAILURE OF THE PROGRAM
TO OPERATE WITH ANY OTHER PROGRAMS), EVEN IF SUCH HOLDER OR OTHER
PARTY HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGES.


