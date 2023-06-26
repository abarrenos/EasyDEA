context('bannerStr')

bannerStrOut <- "\n#     #   #                                     ##           #     #    \n#     #                                          #           #     #    \n#     #  ##             #     #  ####   ## ##    #       #####     #    \n#######   #             #     # #    #   ##  #   #      #    #     #    \n#     #   #             #  #  # #    #   #       #      #    #     #    \n#     #   #             #  #  # #    #   #       #      #    #          \n#     # #####            ## ##   ####   ##      ###      #### #    #    \n                                                                        \n"


test_that('bannerStr() returns a properly formatted string',
    {
	expect_match(bannerStr('Hi world!'), bannerStrOut, all=TRUE)
    }
)

