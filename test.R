#!/usr/bin/env Rscript

print(commandArgs())

# the following worhk only if we are run with 'source()'

fn <- function(x) {
  x + 1 # A comment, kept as part of the source
}                       

# Show the temporary file directory
# where the example was saved

print(getSrcDirectory(fn))
print(getSrcLocation(fn, "line"))

print(getSrcFilename(function(x) {x}, full.names=T))

print(getSrcDirectory(function() {return}))

print(getSrcLocation(function() {}, "line"))


# this exploit Rscript implicit --file argument and should work when
# run as a program and only works if run with Rscript
getScriptPath <- function() {
     
    # this works if we were called with 'source()'
    src.path <- getSrcFilename(function() {}, full.names=T)
    if (! is.null(src.path) && length(src.path) != 0)
            return(src.path)

    # this works for Rscript, it may match more than one path
    # if Rscript was used with several --file= arguments, in
    # which case we will return *only* the first one
    cmd.args <- commandArgs()
    m <- regexpr("(?<=^--file=).+", cmd.args, perl=TRUE)
    script.path <- regmatches(cmd.args, m)
    if (length(script.path) >= 1) return(script.path[1])
    if (length(script.path) == 1) return(script.path)	# may return multiple matches

    # this works for 'R -f', it will return only the first -f argument
    for (i in 1:length(cmd.args) ) {
        print(i); print(cmd.args[i])
        if (cmd.args[i] == '-f') return(cmd.args[i+1])
    }
    
    # if we arrive here, we didn't match anything, turn to last resort
    return (sys.frame(1)$ofile)
}

print("***getScriptPath***")
print(getScriptPath())
print(dirname(getScriptPath()))

mypath <-  function() {
    cmd.args <- commandArgs()
    for ( i in 1:length(cmd.args) ) {
        if (cmd.args[i] == '-f') return(cmd.args[i+1])
    }
    print('last resort')
    return (sys.frame(1)$ofile)
}

print("mypath")
print(mypath())


#
