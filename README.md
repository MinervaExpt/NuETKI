# CC QELike TKI Cross Section for MINERvA
A fork of Andrew's 2021 Cross Section tutorial which serves as the base for my analysis package/framework


Some outstanding framework to dos:
- Rename all of the michel events, michel sidebands stuff to electron events. That's all a holdover from the fact that the tutorial was forked from Mehreen's analysis.

- Rewrite my binning to live in Binning.h and just reference the set of bins I need in runEventLoop, rather than the huge mess I have now

- Rewrite my cuts to allow me to edit the cut value from runEventLoop, it's there for one or two of the cuts but not all of them. Most are hardcoded atm.

- the way I've set up plotting true variables vs reco variables in runEventLoop is very wrong. See if I can figure out how to clean that up...

- More development on grid running... I can submit individual run event loop jobs, but I'd like to set up machinery to submit multiple jobs in parallel and merge the output, that would make running the full FHC data set MUCH easier

- Additional Sideband development, particularly fitting and subtraction. I'm gonna have to look at other analyses for help here. 