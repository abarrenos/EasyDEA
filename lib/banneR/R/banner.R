#' `banneR` mimics UNIX(TM) original banner program in R.
#'
#' `banneR` provides an R `banner()`function that behaves similar to
#' UNIX(TM) `banner` program.
#'
#' @docType package
#'
#' @name banneR
#'
#' @keywords internal
#'
#' @author	José R. Valverde, 2018
"_PACKAGE"

# charScanLine
#
# Convert a character to scan lines
#
# @param ch	character
#
# @param l	the line number of the character
#
# @return	the character scan line corresponding to the specified
#		character and scan line
#
# @author	José R. Valverde
#
charScanLine <- function(ch, l) {
    # This is the representation of each char as a series of 8 8-bit vectors
    #
    # To print a given character we must , for each of the eight lines
    # take the corresponding part and output the resulting horizontal scan 
    # line.
    #
    # This function will return, for a given character, the corresponding
    # horizontal scan line represented as a string of bits, where each
    # hash represents the corresponding bit
    #
    if      (ch == ' ')
        return(c(
        "        ", 
        "        ", 
        "        ", 
        "        ", 
        "        ", 
        "        ", 
        "        ", 
        "        "
	)[l])
    else if (ch == '0')
        return(c(
        " #####  ",
        "#    ## ",
        "#   # # ",
        "#  #  # ",
        "# #   # ",
        "##    # ",
        " #####  ",
        "        "
	)[l])
    else if (ch == "1")
        return(c(
        "     #  ",
        "    ##  ",
        "     #  ",
        "     #  ",
        "     #  ",
        "     #  ",
        "     #  ",
        "        "
        )[l])
    else if (ch == "2")
        return(c(
        " #####  ",
        "#     # ",
        "      # ",
        "  ####  ",
        " #      ",
        "#       ",
        "####### ",
        "        "
        )[l])
    else if (ch == "3")
        return(c(
        " #####  ",
        "#     # ",
        "      # ",
        "  ####  ",
        "      # ",
        "#     # ",
        " #####  ",
        "        "
        )[l])
    else if (ch == "4")
        return(c(
        "#    #  ",
        "#    #  ",
        "#    #  ",
        " ###### ",
        "     #  ",
        "     #  ",
        "     #  ",
        "        "
        )[l])
    else if (ch == "5")
        return(c(
        "####### ",
        "#       ",
        "#       ",
        "######  ",
        "      # ",
        "#     # ",
        " #####  ",
        "        "
        )[l])
    else if (ch == "6")
        return(c(
        " #####  ",
        "#       ",
        "#       ",
        " #####  ",
        "#     # ",
        "#     # ",
        " #####  ",
        "        "
        )[l])
    else if (ch == "7")
        return(c(
        " ###### ",
        "      # ",
        "     #  ",
        "    #   ",
        "   #    ",
        "  #     ",
        " #      ",
        "        "
        )[l])
    else if (ch == "8")
        return(c(
        " #####  ",
        "#     # ",
        "#     # ",
        " #####  ",
        "#     # ",
        "#     # ",
        " #####  ",
        "        "
        )[l])
    else if (ch == "9")
        return(c(
        " #####  ",
        "#     # ",
        "#     # ",
        " #####  ",
        "      # ",
        "      # ",
        " #####  ",
        "        "
        )[l])
    else if (ch == "a")
        return(c(
        "        ",
        "        ",
        " ####   ",
        "     #  ",
        " #####  ",
        "#    #  ",
        " #### # ",
        "        "      
	)[l])
    else if (ch == "b")
        return(c(
        " #      ",
        " #      ",
        " #####  ",
        " #    # ",
        " #    # ",
        " #    # ",
        "# ####  ",
        "        "
        )[l])
    else if (ch == "c")
        return(c(
        "        ",
        "        ",
        "  ###   ",
        " #   #  ",
        "#       ",
        " #   #  ",
        "  ###   ",
        "        "
        )[l])
    else if (ch == "d")
        return(c(
        "     #  ",
        "     #  ",
        " #####  ",
        "#    #  ",
        "#    #  ",
        "#    #  ",
        " #### # ",
        "        "
        )[l])
    else if (ch == "e")
        return(c(
        "        ",
        "        ",
        " ####   ",
        "#    #  ",
        "######  ",
        "#       ",
        " ####   ",
        "        "
        )[l])
    else if (ch == "f")
        return(c(
        "  ###   ",
        " #      ",
        " #      ",
        "####    ",
        " #      ",
        " #      ",
        " #      ",
        "        "        
        )[l])
    else if (ch == "g")
        return(c(
        "        ",
        "        ",
        " #### # ",
        "#    #  ",
        "#    #  ",
        " ####   ",
        "     #  ",
        " ####   "        
        )[l])
    else if (ch == "h")
        return(c(
        " #      ",
        " #      ",
        " ####   ",
        " #   #  ",
        " #   #  ",
        " #   #  ",
        "##   #  ",
        "        "
        )[l])
    else if (ch == "i")
        return(c(
        "  #     ",
        "        ",
        " ##     ",
        "  #     ",
        "  #     ",
        "  #     ",
        "#####   ",
        "        "
        )[l])
    else if (ch == "j")
        return(c(
        "   #    ",
        "        ",
        "  ##    ",
        "   #    ",
        "   #    ",
        "   #    ",
        "#  #    ",
        " ##     "        
        )[l])
    else if (ch == "k")
        return(c(
        " #      ",
        " #  #   ",
        " # #    ",
        " ##     ",
        " # #    ",
        " #  #   ",
        "##  ##  ",
        "        "
        )[l])
    else if (ch == "l")
        return(c(
        "##      ",
        " #      ",
        " #      ",
        " #      ",
        " #      ",
        " #      ",
        "###     ",
        "        "
        )[l])
    else if (ch == "m")
        return(c(
        "        ",
        "        ",
        "# # ##  ",
        "## #  # ",
        "#  #  # ",
        "#  #  # ",
        "#  #  # ",
        "        "
        )[l])
    else if (ch == "n")
        return(c(
        "        ",
        "        ",
        " # ##   ",
        " ##  #  ",
        " #   #  ",
        " #   #  ",
        " #   #  ",
        "        "
        )[l])
    else if (ch == "o")
        return(c(
        "        ",
        "        ",
        " ####   ",
        "#    #  ",
        "#    #  ",
        "#    #  ",
        " ####   ",
        "        "
        )[l])
    else if (ch == "p")
        return(c(
        "        ",
        "        ",
        "#####   ",
        " #   #  ",
        " #   #  ",
        " ####   ",
        " #      ",
        "##      "
        )[l])
    else if (ch == "q")
        return(c(
        "       ",
        "       ",
        " #### #",
        "#    # ",
        "#    # ",
        " ##### ",
        "     # ",
        "     ##" 
        )[l])
    else if (ch == "r")
        return(c(
        "        ",
        "        ",
        "## ##   ",
        " ##  #  ",
        " #      ",
        " #      ",
        "##      ",
        "        "
        )[l])
    else if (ch == "s")
        return(c(
        "        ",
        "        ",
        " ####   ",
        "#       ",
        " ####   ",
        "     #  ",
        "#####   ",
        "        "
        )[l])
    else if (ch == "t")
        return(c(
        "  #     ",
        "  #     ",
        "#####   ",
        "  #     ",
        "  #     ",
        "  #     ",
        "   ##   ",
        "        "
        )[l])
    else if (ch == "u")
        return(c(
        "        ",
        "        ",
        "#    #  ",
        "#    #  ",
        "#    #  ",
        "#    #  ",
        " #### # ",
        "        "        
        )[l])
    else if (ch == "v")
        return(c(
        "        ",
        "        ",
        "#     # ",
        " #   #  ",
        " #   #  ",
        "  # #   ",
        "   #    ",
        "        "        
        )[l])
    else if (ch == "w")
        return(c(
        "        ",
        "        ",
        "#     # ",
        "#     # ",
        "#  #  # ",
        "#  #  # ",
        " ## ##  ",
        "        "               
        )[l])
    else if (ch == "x")
        return(c(
        "        ",
        "        ",
        "#    #  ",
        " #  #   ",
        "  ##    ",
        " #  #   ",
        "#    #  ",
        "        "        
        )[l])
    else if (ch == "y")
        return(c(
        "        ",
        "        ",
        " #   #  ",
        " #   #  ",
        " #   #  ",
        "  ###   ",
        "    #   ",
        " ###    "        
        )[l])
    else if (ch == "z")
        return(c(
        "        ",
        "        ",
        "######  ",
        "    #   ",
        "   #    ",
        "  #     ",
        " ###### ",
        "        "        
        )[l])
    else if (ch == "A")
        return(c(
        "   #    ",
        "  # #   ",
        " #   #  ",
        "#     # ",
        "####### ",
        "#     # ",
        "#     # ",
        "        "
        )[l])
    else if (ch == "B")
        return(c(
        "#####   ",
        "#    #  ",
        "#    #  ",
        "#####   ",
        "#    #  ",
        "#    #  ",
        "#####   ",
        "        "
        )[l])
    else if (ch == "C")
        return(c(
        "  ###   ",
        " #   #  ",
        "#       ",
        "#       ",
        "#       ",
        " #   #  ",
        "  ###   ",
        "        "
        )[l])
    else if (ch == "D")
        return(c(
        "####    ",
        "#   #   ",
        "#    #  ",
        "#    #  ",
        "#    #  ",
        "#   #   ",
        "####    ",
        "        "
        )[l])
    else if (ch == "E")
        return(c(
        "######  ",
        "#       ",
        "#       ",
        "#####   ",
        "#       ",
        "#       ",
        "######  ",
        "        "
        )[l])
    else if (ch == "F")
        return(c(
        "######  ",
        "#       ",
        "#       ",
        "#####   ",
        "#       ",
        "#       ",
        "#       ",
        "        "
        )[l])
    else if (ch == "G")
        return(c(
        "  ###   ",
        " #   #  ",
        "#       ",
        "#       ",
        "#   ### ",
        " #   ## ",
        "  ### # ",
        "        "
        )[l])
    else if (ch == "H")
        return(c(
        "#     # ",
        "#     # ",
        "#     # ",
        "####### ",
        "#     # ",
        "#     # ",
        "#     # ",
        "        "
        )[l])
    else if (ch == "I")
        return(c(
        "#####   ",
        "  #     ",
        "  #     ",
        "  #     ",
        "  #     ",
        "  #     ",
        "#####   ",
        "        "
        )[l])
    else if (ch == "J")
        return(c(
        " #####  ",
        "   #    ",
        "   #    ",
        "   #    ",
        "   #    ",
        " # #    ",
        "  #     ",
        "        "
        )[l])
    else if (ch == "K")
        return(c(
        "#    #  ",
        "#   #   ",
        "#  #    ",
        "###     ",
        "#  #    ",
        "#   #   ",
        "#    #  ",
        "        "        
        )[l])
    else if (ch == "L")
        return(c(
        "#       ",
        "#       ",
        "#       ",
        "#       ",
        "#       ",
        "#       ",
        "######  ",
        "        "        
        )[l])
    else if (ch == "M")
        return(c(
        "##   ## ",
        "# # # # ",
        "#  #  # ",
        "#     # ",
        "#     # ",
        "#     # ",
        "#     # ",
        "        "        
        )[l])
    else if (ch == "N")
        return(c(
        "#     # ",
        "##    # ",
        "# #   # ",
        "#  #  # ",
        "#   # # ",
        "#    ## ",
        "#     # ",
        "        "        
        )[l])
    else if (ch == "O")
        return(c(
        "  ###   ",
        " #   #  ",
        "#     # ",
        "#     # ",
        "#     # ",
        " #   #  ",
        "  ###   ",
        "        "        
        )[l])
    else if (ch == "P")
        return(c(
        "#####   ",
        "#    #  ",
        "#    #  ",
        "#####   ",
        "#       ",
        "#       ",
        "#       ",
        "        "        
        )[l])
    else if (ch == "Q")
        return(c(
        "  ###   ",
        " #   #  ",
        "#     # ",
        "#     # ",
        "# #   # ",
        " # # #  ",
        "  ###   ",
        "     ## "        
        )[l])
    else if (ch == "R")
        return(c(
        "####    ",
        "#   #   ",
        "#   #   ",
        "####    ",
        "#   #   ",
        "#    #  ",
        "#    #  ",
        "        "        
        )[l])
    else if (ch == "S")
        return(c(
        "  ###   ",
        " #   #  ",
        "  #     ",
        "   #    ",
        "    #   ",
        " #   #  ",
        "  ###   ",
        "        "        
        )[l])
    else if (ch == "T")
        return(c(
        "####### ",
        "   #    ",
        "   #    ",
        "   #    ",
        "   #    ",
        "   #    ",
        "   #    ",
        "        "        
        )[l])
    else if (ch == "U")
        return(c(
        "#     # ",
        "#     # ",
        "#     # ",
        "#     # ",
        "#     # ",
        " #   #  ",
        "  ###   ",
        "        "
        )[l])
    else if (ch == "V")
        return(c(
        "#     # ",
        "#     # ",
        " #   #  ",
        " #   #  ",
        "  # #   ",
        "  # #   ",
        "   #    ",
        "        "
        )[l])
    else if (ch == "W")
        return(c(
        "#     # ",
        "#     # ",
        "#     # ",
        "#     # ",
        "#  #  # ",
        "# # # # ",
        " #   #  ",
        "        "
        )[l])
    else if (ch == "X")
        return(c(
        "#     # ",
        "#     # ",
        " #   #  ",
        "  ###   ",
        " #   #  ",
        "#     # ",
        "#     # ",
        "        "
        )[l])
    else if (ch == "Y")
        return(c(
        "#     # ",
        "#     # ",
        " #   #  ",
        "  # #   ",
        "   #    ",
        "   #    ",
        "   #    ",
        "        "
        )[l])
    else if (ch == "Z")
        return(c(
        "####### ",
        "     #  ",
        "    #   ",
        "   #    ",
        "  #     ",
        " #      ",
        "####### ",
        "        "
        )[l])
    else if (ch == ".")
        return(c(
        "        ",
        "        ",
        "        ",
        "        ",
        "        ",
        "        ",
        "   #    ",
        "        "        
        )[l])
    else if (ch == ",")
        return(c(
        "        ",
        "        ",
        "        ",
        "        ",
        "        ",
        "    #   ",
        "    #   ",
        "   #    "        
        )[l])
    else if (ch == ";")
        return(c(
        "        ",
        "        ",
        "        ",
        "    #   ",
        "        ",
        "    #   ",
        "    #   ",
        "   #    "        
        )[l])
    else if (ch == ":")
        return(c(
        "        ",
        "        ",
        "        ",
        "        ",
        "   #    ",
        "        ",
        "   #    ",
        "        "        
        )[l])
    else if (ch == "!")
        return(c(
        "   #    ",
        "   #    ",
        "   #    ",
        "   #    ",
        "   #    ",
        "        ",
        "   #    ",
        "        "        
        )[l])
    else if (ch == "\\")
        return(c(
        "#       ",
        " #      ",
        "  #     ",
        "   #    ",
        "    #   ",
        "     #  ",
        "      # ",
        "        "        
        )[l])
    else if (ch == "/")
        return(c(
        "      # ",
        "     #  ",
        "    #   ",
        "   #    ",
        "  #     ",
        " #      ",
        "#       ",
        "        "        
        )[l])
    else if (ch == "|")
        return(c(
        "   #    ",
        "   #    ",
        "   #    ",
        "   #    ",
        "   #    ",
        "   #    ",
        "   #    ",
        "   #    "        
        )[l])
    else if (ch == "+")
        return(c(
        "        ",
        "   #    ",
        "   #    ",
        " #####  ",
        "   #    ",
        "   #    ",
        "        ",
        "        "        
        )[l])
    else if (ch == "-")
        return(c(
        "        ",
        "        ",
        "        ",
        " #####  ",
        "        ",
        "        ",
        "        ",
        "        "        
        )[l])
    else if (ch == "*")
        return(c(
        "        ",
        " # # #  ",
        "  ###   ",
        "   #    ",
        "  ###   ",
        " # # #  ",
        "        ",
        "        "        
        )[l])
    else if (ch == "=")
        return(c(
        "        ",
        "        ",
        " #####  ",
        "        ",
        " #####  ",
        "        ",
        "        ",
        "        "        
        )[l])
    else if (ch == "_")
        return(c(
	"        ", 
        "        ", 
        "        ", 
        "        ", 
        "        ", 
        "        ", 
        "        ",
        "########"        
        )[l])
    else if (ch == "(")
        return(c(
	"     #  ", 
        "    #   ", 
        "   #    ", 
        "   #    ", 
        "   #    ", 
        "    #   ", 
        "     #  ",
        "        "       
        )[l])
    else if (ch == ")")
        return(c(
	"  #     ", 
        "   #    ", 
        "    #   ", 
        "    #   ", 
        "    #   ", 
        "   #    ", 
        "  #     ",
        "        "       
        )[l])
    else if (ch == "[")
        return(c(
	"   ###  ", 
        "   #    ", 
        "   #    ", 
        "   #    ", 
        "   #    ", 
        "   #    ", 
        "   ###  ",
        "        "       
        )[l])
    else if (ch == "]")
        return(c(
	"  ###   ", 
        "    #   ", 
        "    #   ", 
        "    #   ", 
        "    #   ", 
        "    #   ", 
        "  ###   ",
        "        "       
        )[l])
    else if (ch == "{")
        return(c(
	"    ##  ", 
        "   #    ", 
        "   #    ", 
        "  #     ", 
        "   #    ", 
        "   #    ", 
        "    ##  ",
        "        "       
        )[l])
    else if (ch == "}")
        return(c(
	"  ##    ", 
        "    #   ", 
        "    #   ", 
        "     #  ", 
        "    #   ", 
        "    #   ", 
        "  ##    ",
        "        "       
        )[l])
    else if (ch == "'")
        return(c(
        "    #   ",
        "    #   ",
        "   #    ",       
        "        ",
        "        ",
        "        ",
        "        ",
        "        "
        )[l])
    else if (ch == '"')
        return(c(
        "  # #   ",
        "  # #   ",
        "  # #   ",       
        "        ",
        "        ",
        "        ",
        "        ",
        "        "
        )[l])
    else               
        return(c(
        "########",
        "#      #",
        "#      #",
        "#      #",
        "#      #",
        "#      #",
        "#      #",
        "########"        
        )[l])

}

#' banner
#'
#' \code{banner()} prints out a text in large letters much like UNIX banner
#'	but does not truncate text at 10 characters. Instead, it will print
#'	all the text given in a single "line".
#'
#'	Since output letters are 8 characters wide, and 8 characters high,
#'	this means, in fact, that if you want your text to look good, you 
#'	should truncate it to the ine width limit selected in your R session 
#'	divided by 8.
#'
#'	Maybe in a future release I will modify the package to detect maximum
#'	line width and do the truncation... for now I find it more useful not
#'	to have explicit limits.
#'
#	Banner is used typically to draw attention to special text in your R
#'	output
#'
#' @param text	the text to print (default=" ")
#'
#' @return	nothing
#'
#' @usage	banner(text="something")
#'
#' @examples
#'	banner('_______')
#'	banner()
#'	banner(' Hiya! ')
#'	banner('_______')
#'	\dontshow{ banner('Hello World!') }
#'
#' @author	(C) José R Valverde, CNB-SIC, 2018
#'
#' @license	EU-GPL
#'
#' @export
#
banner <- function(text=" ")
{
    #len <- nchar(text)
    chars <- strsplit(text, "")[[1]]
    
    for (l in 1:8) {
        for (ch in chars) {
            cat(charScanLine(ch, l))
        }
        cat('\n')
    }
}


#' bannerStr
#'
#' \code{bannerStr()} transforms a text for printing in large letters much 
#'	like UNIX banner but does not truncate text at 10 characters. Instead,
#' 	 it will produce all the text given as a single "line".
#'
#	The result of bannerStr() is intended for use with `cat()`
#'
#'	Since output letters are 8 characters wide, and 8 characters high,
#'	this means, in fact, that if you want your text to look good, you 
#'	should truncate it to the ine width limit selected in your R session 
#'	divided by 8.
#'
#'	Maybe in a future release I will modify the package to detect maximum
#'	line width and do the truncation... for now I find it more useful not
#'	to have explicit limits.
#'
#'	Banner is used typically to draw attention to special text in your R
#'	output
#'
#' @param text	the text to print (default=" ")
#'
#' @return	a long character string that spans multiple lines, so that
#'		each line contains the corresponding scan line for the text
#		specified.
#'
#' @usage	cat(bannerStr(text="something"))
#'
#' @examples
#'	cat(bannerStr('Hello world!'))
#'	cat(bannerStr(), sep='')
#'	cat(bannerStr('---------'))
#'	cat(bannerStr('Hi world!'))
#'	cat(bannerStr('---------'))
#'	\dontshow{ cat(banner('Hi World!')) }
#'
#' @author	(C) José R Valverde, CNB-SIC, 2018
#'
#' @license	EU-GPL
#'
#' @export
#
bannerStr <- function(text=" ")
{
    #len <- nchar(text)
    chars <- strsplit(text, "")[[1]]
    
    str=NULL
    for (l in 1:8) {
        str <- paste(str, '\n', sep='')
        for (ch in chars) {
            str <- paste(str, charScanLine(ch, l), sep='')
        }
    }
    str <- paste(str, '\n', sep='')
    return(str)
}



#banner('_______')
#banner()
#banner(' Hiya! ')
#banner('_______')
#banner()
#cat(bannerStr('---------'))
#cat(bannerStr('Hi world!'))
#cat(bannerStr('---------'))
