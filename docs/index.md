--- 
title: "The MumfordBrainStats Mixed Models Series: Companion for the YouTube series"
author: "Jeanette Mumford"
date: "2020-01-31"
site: bookdown::bookdown_site
output: bookdown::gitbook
documentclass: book
cover-image: "./images/mixed_model_book_cover.png"
#bibliography: [book.bib, packages.bib]
biblio-style: apalike
link-citations: yes
#github-repo: rstudio/bookdown-demo
description: ""
---


# Introduction {-}


\includegraphics[width=250px]{./images/mixed_model_book_cover} 

This is a collection of materials that accompanies a [YouTube series on the MumfordBrainStats channel about mixed models](https://www.youtube.com/watch?v=IGHm1XHFWMc&list=PLB2iAtgpI4YEAUiEQ1ZnfMXY-yewNzn9z).  Although I normall focus on material related to neuroimaging, this is for a general audience.  Each of these chapters should be understandable without watching the video, but one would probably gain the most by watching the videos as well.  The chapter titles indicate which video in the series goes along with that chapter.   Not all videos have chapter (yet), since I'm only including chapters with code for now.  In the future I may try to smooth this out more to read more like an actual book.

Although this series isn't necessarily set up to be a tutorial on how to run a mixed model properly, it would still be useful for somebody who is new to mixed models.  I'm focusing on odds and ends that I've seen come up that I feel are less frequently discussed.  Basically mistakes I've made or seen that I'd like to help you avoid.  Interestingly, many of these mistakes are made by those with some knowledge of mixed models.  I think they often arises from  a "knows enough to be dangerous" situation.  Do not misinterpret that as it is good to be dangerous, brave, make mistakes and then fix them.  That's how we learn!

Specific goals of this series:

* Make it easier for you to conceptualize what a mixed model is doing
* Describe the two stage random effects formulation and corresponding two stage summary statistics approach for modeling repeated measures data
* Clarify when regularization is occurring and what it is doing
* Clarifying what a conditional mode (BLUP) is and how they should *not* be used
* When a two stage summary statistics approach is likely okay to use

This series is still in development.  Future topics will involve: crossed/nested random effects, tips for planning out your analysis setup, simplifying random effects structures (when is it okay) and hopefully some examples from the audience of this series about their trials and tribulations with mixed models.  I'm also looking for solid ways to test model assumptions.  So far I haven't found anything super satisfying.

For feedback, you can find me on twitter: @mumbrainstats.
