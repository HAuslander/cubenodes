# the following quick and dirty program
# acts on output from the cubenodes programs that has
# most of the lines deleted as so (so that remaining lines start with "d" or "H"):

# grep -v Div blah32.txt > blah31.txt
# grep -v Well blah31.txt > blah32.txt
# 
# e.g.
"""


NDIV 0 280
Hist update for nnode count 624
d 1.0 3024
d 3.0 576
NDIV 1 164
Hist update for nnode count 680
d 1.0 3000
d 3.0 864
NDIV 2 492
Hist update for nnode count 736
d 1.0 2976
d 3.0 1152
NDIV 3 324
Hist update for nnode count 792
d 1.0 2952
d 3.0 1440
NDIV 4 230
Hist update for nnode count 846
d 1.0 2930
d 3.0 1722
NDIV 5 545
Hist update for nnode count 902
d 1.0 2930
d 3.0 1698
d 9.0 288
NDIV 6 32
Hist update for n


"""

# the initial script being run for the following is: %run  /Users/hrvojehrgovcic/quant/latticegas_cubenodes14.py  -t 100000  --xprob 1.0 --maxnode 80000 --dim 3


from copy import copy


import matplotlib.pyplot as plt
import random as rn
#from bisect import bisect_left
f = open("blah32.txt")
flines = f.readlines()
f.close()
x = []
d1 = []
d3 = []
d9 = []
d27 = []

thisx = 0
this1 = 0
this3 = 0
this9 = 0
this27 = 0
for il in flines:
    if il[0] == "H":
        x.append(thisx)
        d1.append(this1)
        d3.append(this3)
        d9.append(this9)
        d27.append(this27)
        elem = il.strip().split(' ')
        thisx =  int(elem[-1])
    elif il.find("d 1.0") == 0:
        elem = il.strip().split(' ')
        this1 = int(elem[2])

    elif il.find("d 3.0") == 0:
        elem = il.strip().split(' ')
        this3 = int(elem[2])

    elif il.find("d 9.0") == 0:
        elem = il.strip().split(' ')
        this9 = int(elem[2])

    elif il.find("d 27.0") == 0:
        elem = il.strip().split(' ')
        this27 = int(elem[2])

# plot lines
plt.plot(x, d1, label = "d1")
plt.plot(x, d3, label = "d3")
plt.plot(x, d9, label = "d9")
plt.plot(x, d27, label = "d27")
plt.legend()
plt.show()

"""


0
NDIV 0 280
Hist update for nnode count 624
d 1.0 3024
d 3.0 576
NDIV 1 164
Hist update for nnode count 680
d 1.0 3000
d 3.0 864
NDIV 2 492
Hist update for nnode count 736
d 1.0 2976
d 3.0 1152
NDIV 3 324
Hist update for nnode count 792
d 1.0 2952
d 3.0 1440
NDIV 4 230
Hist update for nnode count 846
d 1.0 2930
d 3.0 1722
NDIV 5 545
Hist update for nnode count 902
d 1.0 2930
d 3.0 1698
d 9.0 288
NDIV 6 32
Hist update for nnode count 958
d 1.0 2906
d 3.0 1986
d 9.0 288
NDIV 7 367
Hist update for nnode count 1012
d 1.0 2884
d 3.0 2268
d 9.0 288
NDIV 8 197
Hist update for nnode count 1068
d 1.0 2860
d 3.0 2556
d 9.0 288
NDIV 9 191
Hist update for nnode count 1122
d 1.0 2838
d 3.0 2838
d 9.0 288
NDIV 10 366
Hist update for nnode count 1172
d 1.0 2820
d 3.0 3108
d 9.0 288
NDIV 11 234
Hist update for nnode count 1226
d 1.0 2798
d 3.0 3390
d 9.0 288
NDIV 12 182
Hist update for nnode count 1270
d 1.0 2782
d 3.0 3630
d 9.0 288
NDIV 13 125
Hist update for nnode count 1324
d 1.0 2760
d 3.0 3912
d 9.0 288
NDIV 14 427
Hist update for nnode count 1378
d 1.0 2738
d 3.0 4194
d 9.0 288
NDIV 15 97
Hist update for nnode count 1432
d 1.0 2716
d 3.0 4476
d 9.0 288
NDIV 16 355
Hist update for nnode count 1486
d 1.0 2694
d 3.0 4758
d 9.0 288
NDIV 17 430
Hist update for nnode count 1540
d 1.0 2672
d 3.0 5040
d 9.0 288
NDIV 18 86
Hist update for nnode count 1594
d 1.0 2650
d 3.0 5322
d 9.0 288
NDIV 19 364
Hist update for nnode count 1648
d 1.0 2628
d 3.0 5604
d 9.0 288
NDIV 20 17
Hist update for nnode count 1702
d 1.0 2606
d 3.0 5886
d 9.0 288
NDIV 21 453
Hist update for nnode count 1756
d 1.0 2584
d 3.0 6168
d 9.0 288
NDIV 22 28
Hist update for nnode count 1812
d 1.0 2560
d 3.0 6456
d 9.0 288
NDIV 23 357
Hist update for nnode count 1864
d 1.0 2540
d 3.0 6732
d 9.0 288
NDIV 24 132
Hist update for nnode count 1908
d 1.0 2524
d 3.0 6972
d 9.0 288
NDIV 25 217
Hist update for nnode count 1960
d 1.0 2504
d 3.0 7248
d 9.0 288
NDIV 26 117
Hist update for nnode count 2004
d 1.0 2488
d 3.0 7488
d 9.0 288
NDIV 27 423
Hist update for nnode count 2050
d 1.0 2474
d 3.0 7746
d 9.0 288
NDIV 28 250
Hist update for nnode count 2106
d 1.0 2450
d 3.0 8034
d 9.0 288
NDIV 29 101
Hist update for nnode count 2160
d 1.0 2428
d 3.0 8316
d 9.0 288
NDIV 30 398
Hist update for nnode count 2204
d 1.0 2412
d 3.0 8556
d 9.0 288
NDIV 31 67
Hist update for nnode count 2258
d 1.0 2390
d 3.0 8838
d 9.0 288
NDIV 32 472
Hist update for nnode count 2312
d 1.0 2368
d 3.0 9120
d 9.0 288
NDIV 33 298
Hist update for nnode count 2366
d 1.0 2346
d 3.0 9402
d 9.0 288
NDIV 34 350
Hist update for nnode count 2400
d 1.0 2336
d 3.0 9600
d 9.0 288
NDIV 35 370
Hist update for nnode count 2444
d 1.0 2320
d 3.0 9840
d 9.0 288
NDIV 36 114
Hist update for nnode count 2498
d 1.0 2298
d 3.0 10122
d 9.0 288
NDIV 37 58
Hist update for nnode count 2552
d 1.0 2276
d 3.0 10404
d 9.0 288
NDIV 38 187
Hist update for nnode count 2596
d 1.0 2260
d 3.0 10644
d 9.0 288
NDIV 39 388
Hist update for nnode count 2630
d 1.0 2250
d 3.0 10842
d 9.0 288
NDIV 40 360
Hist update for nnode count 2674
d 1.0 2234
d 3.0 11082
d 9.0 288
NDIV 41 396
Hist update for nnode count 2708
d 1.0 2224
d 3.0 11280
d 9.0 288
NDIV 42 199
Hist update for nnode count 2752
d 1.0 2208
d 3.0 11520
d 9.0 288
NDIV 43 303
Hist update for nnode count 2786
d 1.0 2198
d 3.0 11718
d 9.0 288
NDIV 44 219
Hist update for nnode count 2812
d 1.0 2192
d 3.0 11880
d 9.0 288
NDIV 45 285
Hist update for nnode count 2846
d 1.0 2182
d 3.0 12078
d 9.0 288
NDIV 46 371
Hist update for nnode count 2888
d 1.0 2168
d 3.0 12312
d 9.0 288
NDIV 47 483
Hist update for nnode count 2940
d 1.0 2148
d 3.0 12588
d 9.0 288
NDIV 48 68
Hist update for nnode count 2982
d 1.0 2134
d 3.0 12822
d 9.0 288
NDIV 49 74
Hist update for nnode count 3026
d 1.0 2118
d 3.0 13062
d 9.0 288
NDIV 50 200
Hist update for nnode count 3080
d 1.0 2096
d 3.0 13344
d 9.0 288
NDIV 51 21
Hist update for nnode count 3132
d 1.0 2076
d 3.0 13620
d 9.0 288
NDIV 52 145
Hist update for nnode count 3184
d 1.0 2056
d 3.0 13896
d 9.0 288
NDIV 53 411
Hist update for nnode count 3236
d 1.0 2036
d 3.0 14172
d 9.0 288
NDIV 54 246
Hist update for nnode count 3290
d 1.0 2014
d 3.0 14454
d 9.0 288
NDIV 55 267
Hist update for nnode count 3342
d 1.0 1994
d 3.0 14730
d 9.0 288
NDIV 56 256
Hist update for nnode count 3394
d 1.0 1974
d 3.0 15006
d 9.0 288
NDIV 57 470
Hist update for nnode count 3444
d 1.0 1956
d 3.0 15276
d 9.0 288
NDIV 58 16
Hist update for nnode count 3488
d 1.0 1940
d 3.0 15516
d 9.0 288
NDIV 59 102
Hist update for nnode count 3530
d 1.0 1926
d 3.0 15750
d 9.0 288
NDIV 60 76
Hist update for nnode count 3574
d 1.0 1910
d 3.0 15990
d 9.0 288
NDIV 61 54
Hist update for nnode count 3626
d 1.0 1890
d 3.0 16266
d 9.0 288
NDIV 62 94
Hist update for nnode count 3670
d 1.0 1874
d 3.0 16506
d 9.0 288
NDIV 63 341
Hist update for nnode count 3722
d 1.0 1854
d 3.0 16782
d 9.0 288
NDIV 64 362
Hist update for nnode count 3754
d 1.0 1846
d 3.0 16974
d 9.0 288
NDIV 65 495
Hist update for nnode count 3794
d 1.0 1834
d 3.0 17202
d 9.0 288
NDIV 66 409
Hist update for nnode count 3846
d 1.0 1814
d 3.0 17478
d 9.0 288
NDIV 67 425
Hist update for nnode count 3880
d 1.0 1804
d 3.0 17676
d 9.0 288
NDIV 68 451
Hist update for nnode count 3922
d 1.0 1790
d 3.0 17910
d 9.0 288
NDIV 69 442
Hist update for nnode count 3974
d 1.0 1770
d 3.0 18186
d 9.0 288
NDIV 70 447
Hist update for nnode count 4018
d 1.0 1754
d 3.0 18426
d 9.0 288
NDIV 71 93
Hist update for nnode count 4042
d 1.0 1750
d 3.0 18582
d 9.0 288
NDIV 72 401
Hist update for nnode count 4096
d 1.0 1728
d 3.0 18864
d 9.0 288
NDIV 73 274
Hist update for nnode count 4138
d 1.0 1714
d 3.0 19098
d 9.0 288
NDIV 74 507
Hist update for nnode count 4192
d 1.0 1692
d 3.0 19380
d 9.0 288
NDIV 75 422
Hist update for nnode count 4230
d 1.0 1682
d 3.0 19602
d 9.0 288
NDIV 76 429
Hist update for nnode count 4274
d 1.0 1670
d 3.0 19854
d 9.0 288
NDIV 77 152
Hist update for nnode count 4328
d 1.0 1648
d 3.0 20136
d 9.0 288
NDIV 78 198
Hist update for nnode count 4360
d 1.0 1640
d 3.0 20328
d 9.0 288
NDIV 79 320
Hist update for nnode count 4412
d 1.0 1620
d 3.0 20604
d 9.0 288
NDIV 80 392
Hist update for nnode count 4464
d 1.0 1600
d 3.0 20880
d 9.0 288
NDIV 81 180
Hist update for nnode count 4516
d 1.0 1580
d 3.0 21156
d 9.0 288
NDIV 82 60
Hist update for nnode count 4554
d 1.0 1570
d 3.0 21378
d 9.0 288
NDIV 83 195
Hist update for nnode count 4602
d 1.0 1554
d 3.0 21642
d 9.0 288
NDIV 84 34
Hist update for nnode count 4654
d 1.0 1534
d 3.0 21918
d 9.0 288
NDIV 85 116
Hist update for nnode count 4692
d 1.0 1524
d 3.0 22140
d 9.0 288
NDIV 86 460
Hist update for nnode count 4736
d 1.0 1508
d 3.0 22380
d 9.0 288
NDIV 87 456
Hist update for nnode count 4790
d 1.0 1486
d 3.0 22662
d 9.0 288
NDIV 88 501
Hist update for nnode count 4836
d 1.0 1472
d 3.0 22920
d 9.0 288
NDIV 89 477
Hist update for nnode count 4852
d 1.0 1472
d 3.0 23040
d 9.0 288
NDIV 90 353
Hist update for nnode count 4882
d 1.0 1466
d 3.0 23226
d 9.0 288
NDIV 91 289
Hist update for nnode count 4922
d 1.0 1454
d 3.0 23454
d 9.0 288
NDIV 92 420
Hist update for nnode count 4958
d 1.0 1446
d 3.0 23670
d 9.0 288
NDIV 93 10
Hist update for nnode count 5000
d 1.0 1432
d 3.0 23904
d 9.0 288
NDIV 94 266
Hist update for nnode count 5032
d 1.0 1424
d 3.0 24096
d 9.0 288
NDIV 95 254
Hist update for nnode count 5064
d 1.0 1416
d 3.0 24288
d 9.0 288
NDIV 96 161
Hist update for nnode count 5098
d 1.0 1406
d 3.0 24486
d 9.0 288
NDIV 97 476
Hist update for nnode count 5144
d 1.0 1392
d 3.0 24744
d 9.0 288
NDIV 98 231
Hist update for nnode count 5160
d 1.0 1392
d 3.0 24864
d 9.0 288
NDIV 99 70
Hist update for nnode count 5204
d 1.0 1376
d 3.0 25104
d 9.0 288
NDIV 100 193
Hist update for nnode count 5256
d 1.0 1356
d 3.0 25380
d 9.0 288
NDIV 101 129
Hist update for nnode count 5300
d 1.0 1340
d 3.0 25620
d 9.0 288
NDIV 102 502
Hist update for nnode count 5344
d 1.0 1328
d 3.0 25872
d 9.0 288
NDIV 103 317
Hist update for nnode count 5398
d 1.0 1306
d 3.0 26154
d 9.0 288
NDIV 104 5
Hist update for nnode count 5440
d 1.0 1292
d 3.0 26388
d 9.0 288
NDIV 105 318
Hist update for nnode count 5478
d 1.0 1282
d 3.0 26610
d 9.0 288
NDIV 106 343
Hist update for nnode count 5526
d 1.0 1266
d 3.0 26874
d 9.0 288
NDIV 107 167
Hist update for nnode count 5558
d 1.0 1258
d 3.0 27066
d 9.0 288
NDIV 108 14
Hist update for nnode count 5598
d 1.0 1246
d 3.0 27294
d 9.0 288
NDIV 109 185
Hist update for nnode count 5642
d 1.0 1230
d 3.0 27534
d 9.0 288
NDIV 110 135
Hist update for nnode count 5672
d 1.0 1224
d 3.0 27720
d 9.0 288
NDIV 111 92
Hist update for nnode count 5698
d 1.0 1218
d 3.0 27882
d 9.0 288
NDIV 112 103
Hist update for nnode count 5742
d 1.0 1202
d 3.0 28122
d 9.0 288
NDIV 113 64
Hist update for nnode count 5784
d 1.0 1188
d 3.0 28356
d 9.0 288
NDIV 114 47
Hist update for nnode count 5838
d 1.0 1166
d 3.0 28638
d 9.0 288
NDIV 115 40
Hist update for nnode count 5888
d 1.0 1148
d 3.0 28908
d 9.0 288
NDIV 116 225
Hist update for nnode count 5926
d 1.0 1138
d 3.0 29130
d 9.0 288
NDIV 117 295
Hist update for nnode count 5952
d 1.0 1132
d 3.0 29292
d 9.0 288
Found alternative choice -- moving from  550 to 210
NDIV 118 550
Hist update for nnode count 5984
d 1.0 1124
d 3.0 29484
d 9.0 288
NDIV 119 130
Hist update for nnode count 6014
d 1.0 1118
d 3.0 29670
d 9.0 288
NDIV 120 51
Hist update for nnode count 6040
d 1.0 1112
d 3.0 29832
d 9.0 288
NDIV 121 434
Hist update for nnode count 6082
d 1.0 1098
d 3.0 30066
d 9.0 288
NDIV 122 323
Hist update for nnode count 6116
d 1.0 1088
d 3.0 30264
d 9.0 288
Found alternative choice -- moving from  549 to 209
NDIV 123 549
Hist update for nnode count 6150
d 1.0 1078
d 3.0 30462
d 9.0 288
NDIV 124 276
Hist update for nnode count 6204
d 1.0 1056
d 3.0 30744
d 9.0 288
NDIV 125 278
Hist update for nnode count 6256
d 1.0 1036
d 3.0 31020
d 9.0 288
NDIV 126 440
Hist update for nnode count 6306
d 1.0 1018
d 3.0 31290
d 9.0 288
NDIV 127 1
Hist update for nnode count 6334
d 1.0 1014
d 3.0 31470
d 9.0 288
NDIV 128 269
Hist update for nnode count 6358
d 1.0 1010
d 3.0 31626
d 9.0 288
NDIV 129 263
Hist update for nnode count 6390
d 1.0 1002
d 3.0 31818
d 9.0 288
NDIV 130 342
Hist update for nnode count 6424
d 1.0 996
d 3.0 32028
d 9.0 288
NDIV 131 414
Hist update for nnode count 6448
d 1.0 992
d 3.0 32184
d 9.0 288
NDIV 132 33
Hist update for nnode count 6488
d 1.0 980
d 3.0 32412
d 9.0 288
NDIV 133 158
Hist update for nnode count 6538
d 1.0 962
d 3.0 32682
d 9.0 288
NDIV 134 43
Hist update for nnode count 6576
d 1.0 952
d 3.0 32904
d 9.0 288
NDIV 135 407
Hist update for nnode count 6618
d 1.0 938
d 3.0 33138
d 9.0 288
NDIV 136 384
Hist update for nnode count 6660
d 1.0 924
d 3.0 33372
d 9.0 288
NDIV 137 404
Hist update for nnode count 6700
d 1.0 912
d 3.0 33600
d 9.0 288
NDIV 138 431
Hist update for nnode count 6728
d 1.0 908
d 3.0 33780
d 9.0 288
NDIV 139 271
Hist update for nnode count 6746
d 1.0 906
d 3.0 33906
d 9.0 288
NDIV 140 497
Hist update for nnode count 6788
d 1.0 892
d 3.0 34140
d 9.0 288
NDIV 141 208
Hist update for nnode count 6842
d 1.0 870
d 3.0 34422
d 9.0 288
Found alternative choice -- moving from  546 to 138
NDIV 142 546
Hist update for nnode count 6892
d 1.0 852
d 3.0 34692
d 9.0 288
NDIV 143 322
Hist update for nnode count 6916
d 1.0 848
d 3.0 34848
d 9.0 288
NDIV 144 249
Hist update for nnode count 6958
d 1.0 834
d 3.0 35082
d 9.0 288
NDIV 145 310
Hist update for nnode count 6998
d 1.0 822
d 3.0 35310
d 9.0 288
Found alternative choice -- moving from  537 to 6129
NDIV 146 537
Hist update for nnode count 7038
d 1.0 810
d 3.0 35538
d 9.0 288
NDIV 147 171
Hist update for nnode count 7072
d 1.0 800
d 3.0 35736
d 9.0 288
NDIV 148 127
Hist update for nnode count 7120
d 1.0 784
d 3.0 36000
d 9.0 288
NDIV 149 163
Hist update for nnode count 7168
d 1.0 768
d 3.0 36264
d 9.0 288
NDIV 150 439
Hist update for nnode count 7196
d 1.0 764
d 3.0 36444
d 9.0 288
NDIV 151 207
Hist update for nnode count 7228
d 1.0 756
d 3.0 36636
d 9.0 288
NDIV 152 290
Hist update for nnode count 7262
d 1.0 750
d 3.0 36846
d 9.0 288
Found alternative choice -- moving from  529 to 541
NDIV 153 529
Hist update for nnode count 7294
d 1.0 742
d 3.0 37038
d 9.0 288
NDIV 154 122
Hist update for nnode count 7326
d 1.0 738
d 3.0 37242
d 9.0 288
NDIV 155 351
Hist update for nnode count 7334
d 1.0 738
d 3.0 37314
d 9.0 288
NDIV 156 99
Hist update for nnode count 7370
d 1.0 730
d 3.0 37530
d 9.0 288
NDIV 157 142
Hist update for nnode count 7402
d 1.0 722
d 3.0 37722
d 9.0 288
NDIV 158 0
Hist update for nnode count 7426
d 1.0 718
d 3.0 37878
d 9.0 288
NDIV 159 107
Hist update for nnode count 7444
d 1.0 716
d 3.0 38004
d 9.0 288
NDIV 160 242
Hist update for nnode count 7478
d 1.0 706
d 3.0 38202
d 9.0 288
NDIV 161 406
Hist update for nnode count 7494
d 1.0 706
d 3.0 38322
d 9.0 288
NDIV 162 222
Hist update for nnode count 7518
d 1.0 702
d 3.0 38478
d 9.0 288
NDIV 163 239
Hist update for nnode count 7570
d 1.0 682
d 3.0 38754
d 9.0 288
NDIV 164 216
Hist update for nnode count 7598
d 1.0 678
d 3.0 38934
d 9.0 288
NDIV 165 309
Hist update for nnode count 7616
d 1.0 676
d 3.0 39060
d 9.0 288
NDIV 166 245
Hist update for nnode count 7652
d 1.0 668
d 3.0 39276
d 9.0 288
NDIV 167 471
Hist update for nnode count 7688
d 1.0 660
d 3.0 39492
d 9.0 288
NDIV 168 72
Hist update for nnode count 7722
d 1.0 650
d 3.0 39690
d 9.0 288
NDIV 169 45
Hist update for nnode count 7746
d 1.0 646
d 3.0 39846
d 9.0 288
NDIV 170 315
Hist update for nnode count 7784
d 1.0 636
d 3.0 40068
d 9.0 288
NDIV 171 284
Hist update for nnode count 7826
d 1.0 622
d 3.0 40302
d 9.0 288
NDIV 172 226
Hist update for nnode count 7856
d 1.0 616
d 3.0 40488
d 9.0 288
NDIV 173 349
Hist update for nnode count 7882
d 1.0 610
d 3.0 40650
d 9.0 288
NDIV 174 293
Hist update for nnode count 7930
d 1.0 594
d 3.0 40914
d 9.0 288
NDIV 175 6
Hist update for nnode count 7968
d 1.0 584
d 3.0 41136
d 9.0 288
NDIV 176 302
Hist update for nnode count 8006
d 1.0 574
d 3.0 41358
d 9.0 288
NDIV 177 38
Hist update for nnode count 8018
d 1.0 574
d 3.0 41454
d 9.0 288
NDIV 178 381
Hist update for nnode count 8066
d 1.0 558
d 3.0 41718
d 9.0 288
NDIV 179 172
Hist update for nnode count 8098
d 1.0 550
d 3.0 41910
d 9.0 288
NDIV 180 59
Hist update for nnode count 8120
d 1.0 548
d 3.0 42060
d 9.0 288
NDIV 181 179
Hist update for nnode count 8146
d 1.0 546
d 3.0 42234
d 9.0 288
NDIV 182 488
Hist update for nnode count 8158
d 1.0 546
d 3.0 42330
d 9.0 288
NDIV 183 332
Hist update for nnode count 8180
d 1.0 544
d 3.0 42480
d 9.0 288
NDIV 184 413
Hist update for nnode count 8198
d 1.0 542
d 3.0 42606
d 9.0 288
NDIV 185 458
Hist update for nnode count 8240
d 1.0 528
d 3.0 42840
d 9.0 288
NDIV 186 372
Hist update for nnode count 8276
d 1.0 520
d 3.0 43056
d 9.0 288
NDIV 187 168
Hist update for nnode count 8316
d 1.0 508
d 3.0 43284
d 9.0 288
NDIV 188 368
Hist update for nnode count 8354
d 1.0 498
d 3.0 43506
d 9.0 288
NDIV 189 178
Hist update for nnode count 8376
d 1.0 496
d 3.0 43656
d 9.0 288
NDIV 190 166
Hist update for nnode count 8392
d 1.0 496
d 3.0 43776
d 9.0 288
NDIV 191 18
Hist update for nnode count 8428
d 1.0 488
d 3.0 43992
d 9.0 288
NDIV 192 184
Hist update for nnode count 8452
d 1.0 484
d 3.0 44148
d 9.0 288
NDIV 193 356
Hist update for nnode count 8468
d 1.0 484
d 3.0 44268
d 9.0 288
NDIV 194 27
Hist update for nnode count 8492
d 1.0 480
d 3.0 44424
d 9.0 288
NDIV 195 224
Hist update for nnode count 8520
d 1.0 476
d 3.0 44604
d 9.0 288
NDIV 196 479
Hist update for nnode count 8544
d 1.0 472
d 3.0 44760
d 9.0 288
NDIV 197 326
Hist update for nnode count 8570
d 1.0 466
d 3.0 44922
d 9.0 288
NDIV 198 63
Hist update for nnode count 8602
d 1.0 458
d 3.0 45114
d 9.0 288
NDIV 199 11
Hist update for nnode count 8626
d 1.0 454
d 3.0 45270
d 9.0 288
NDIV 200 259
Hist update for nnode count 8658
d 1.0 446
d 3.0 45462
d 9.0 288
NDIV 201 511
Hist update for nnode count 8694
d 1.0 438
d 3.0 45678
d 9.0 288
NDIV 202 474
Hist update for nnode count 8736
d 1.0 424
d 3.0 45912
d 9.0 288
NDIV 203 151
Hist update for nnode count 8756
d 1.0 424
d 3.0 46056
d 9.0 288
NDIV 204 238
Hist update for nnode count 8776
d 1.0 424
d 3.0 46200
d 9.0 288
NDIV 205 291
Hist update for nnode count 8806
d 1.0 418
d 3.0 46386
d 9.0 288
Found alternative choice -- moving from  533 to 82
NDIV 206 533
Hist update for nnode count 8826
d 1.0 418
d 3.0 46530
d 9.0 288
NDIV 207 181
Hist update for nnode count 8856
d 1.0 412
d 3.0 46716
d 9.0 288
NDIV 208 390
Hist update for nnode count 8880
d 1.0 408
d 3.0 46872
d 9.0 288
NDIV 209 84
Hist update for nnode count 8920
d 1.0 396
d 3.0 47100
d 9.0 288
NDIV 210 426
Hist update for nnode count 8936
d 1.0 396
d 3.0 47220
d 9.0 288
NDIV 211 124
Hist update for nnode count 8948
d 1.0 396
d 3.0 47316
d 9.0 288
NDIV 212 57
Hist update for nnode count 8986
d 1.0 386
d 3.0 47538
d 9.0 288
NDIV 213 397
Hist update for nnode count 9010
d 1.0 382
d 3.0 47694
d 9.0 288
NDIV 214 465
Hist update for nnode count 9042
d 1.0 374
d 3.0 47886
d 9.0 288
NDIV 215 9
Hist update for nnode count 9060
d 1.0 372
d 3.0 48012
d 9.0 288
NDIV 216 37
Hist update for nnode count 9098
d 1.0 362
d 3.0 48234
d 9.0 288
NDIV 217 75
Hist update for nnode count 9116
d 1.0 360
d 3.0 48360
d 9.0 288
NDIV 218 493
Hist update for nnode count 9132
d 1.0 360
d 3.0 48480
d 9.0 288
NDIV 219 106
Hist update for nnode count 9158
d 1.0 354
d 3.0 48642
d 9.0 288
NDIV 220 405
Hist update for nnode count 9170
d 1.0 354
d 3.0 48738
d 9.0 288
NDIV 221 133
Hist update for nnode count 9190
d 1.0 354
d 3.0 48882
d 9.0 288
NDIV 222 148
Hist update for nnode count 9212
d 1.0 352
d 3.0 49032
d 9.0 288
NDIV 223 131
Hist update for nnode count 9228
d 1.0 352
d 3.0 49152
d 9.0 288
NDIV 224 25
Hist update for nnode count 9252
d 1.0 348
d 3.0 49308
d 9.0 288
NDIV 225 531
Hist update for nnode count 9306
d 1.0 348
d 3.0 49286
d 9.0 570
NDIV 226 336
Hist update for nnode count 9334
d 1.0 344
d 3.0 49466
d 9.0 570
NDIV 227 424
Hist update for nnode count 9360
d 1.0 342
d 3.0 49640
d 9.0 570
NDIV 228 36
Hist update for nnode count 9398
d 1.0 336
d 3.0 49874
d 9.0 570
Found alternative choice -- moving from  547 to 546
NDIV 229 547
Hist update for nnode count 9452
d 1.0 336
d 3.0 49852
d 9.0 852
NDIV 230 402
Hist update for nnode count 9486
d 1.0 326
d 3.0 50050
d 9.0 852
NDIV 231 344
Hist update for nnode count 9508
d 1.0 324
d 3.0 50200
d 9.0 852
NDIV 232 19
Hist update for nnode count 9534
d 1.0 318
d 3.0 50362
d 9.0 852
NDIV 233 354
Hist update for nnode count 9552
d 1.0 316
d 3.0 50488
d 9.0 852
NDIV 234 24
Hist update for nnode count 9576
d 1.0 312
d 3.0 50644
d 9.0 852
NDIV 235 347
Hist update for nnode count 9608
d 1.0 304
d 3.0 50836
d 9.0 852
NDIV 236 186
Hist update for nnode count 9624
d 1.0 304
d 3.0 50956
d 9.0 852
NDIV 237 340
Hist update for nnode count 9644
d 1.0 304
d 3.0 51100
d 9.0 852
NDIV 238 441
Hist update for nnode count 9684
d 1.0 292
d 3.0 51328
d 9.0 852
NDIV 239 194
Hist update for nnode count 9732
d 1.0 276
d 3.0 51592
d 9.0 852
NDIV 240 334
Hist update for nnode count 9760
d 1.0 272
d 3.0 51772
d 9.0 852
NDIV 241 153
Hist update for nnode count 9780
d 1.0 272
d 3.0 51916
d 9.0 852
NDIV 242 491
Hist update for nnode count 9798
d 1.0 270
d 3.0 52042
d 9.0 852
NDIV 243 454
Hist update for nnode count 9842
d 1.0 258
d 3.0 52294
d 9.0 852
NDIV 244 90
Hist update for nnode count 9854
d 1.0 258
d 3.0 52390
d 9.0 852
NDIV 245 189
Hist update for nnode count 9862
d 1.0 258
d 3.0 52462
d 9.0 852
NDIV 246 95
Hist update for nnode count 9878
d 1.0 258
d 3.0 52582
d 9.0 852
NDIV 247 331
Hist update for nnode count 9894
d 1.0 258
d 3.0 52702
d 9.0 852
NDIV 248 240
Hist update for nnode count 9916
d 1.0 256
d 3.0 52852
d 9.0 852
NDIV 249 190
Hist update for nnode count 9932
d 1.0 256
d 3.0 52972
d 9.0 852
NDIV 250 121
Hist update for nnode count 9948
d 1.0 256
d 3.0 53092
d 9.0 852
NDIV 251 98
Hist update for nnode count 9964
d 1.0 256
d 3.0 53212
d 9.0 852
NDIV 252 218
Hist update for nnode count 9986
d 1.0 254
d 3.0 53362
d 9.0 852
NDIV 253 450
Hist update for nnode count 10004
d 1.0 252
d 3.0 53488
d 9.0 852
NDIV 254 111
Hist update for nnode count 10046
d 1.0 242
d 3.0 53734
d 9.0 852
NDIV 255 206
Hist update for nnode count 10066
d 1.0 242
d 3.0 53878
d 9.0 852
NDIV 256 110
Hist update for nnode count 10084
d 1.0 240
d 3.0 54004
d 9.0 852
NDIV 257 475
Hist update for nnode count 10104
d 1.0 240
d 3.0 54148
d 9.0 852
NDIV 258 89
Hist update for nnode count 10120
d 1.0 240
d 3.0 54268
d 9.0 852
NDIV 259 466
Hist update for nnode count 10136
d 1.0 240
d 3.0 54388
d 9.0 852
NDIV 260 140
Hist update for nnode count 10158
d 1.0 238
d 3.0 54538
d 9.0 852
NDIV 261 438
Hist update for nnode count 10184
d 1.0 232
d 3.0 54700
d 9.0 852
NDIV 262 159
Hist update for nnode count 10208
d 1.0 228
d 3.0 54856
d 9.0 852
NDIV 263 157
Hist update for nnode count 10242
d 1.0 218
d 3.0 55054
d 9.0 852
NDIV 264 44
Hist update for nnode count 10258
d 1.0 218
d 3.0 55174
d 9.0 852
NDIV 265 301
Hist update for nnode count 10288
d 1.0 212
d 3.0 55360
d 9.0 852
NDIV 266 134
Hist update for nnode count 10306
d 1.0 210
d 3.0 55486
d 9.0 852
NDIV 267 335
Hist update for nnode count 10318
d 1.0 210
d 3.0 55582
d 9.0 852
NDIV 268 176
Hist update for nnode count 10342
d 1.0 206
d 3.0 55738
d 9.0 852
NDIV 269 526
Hist update for nnode count 10396
d 1.0 206
d 3.0 55716
d 9.0 1134
NDIV 270 437
Hist update for nnode count 10420
d 1.0 202
d 3.0 55872
d 9.0 1134
NDIV 271 258
Hist update for nnode count 10444
d 1.0 198
d 3.0 56028
d 9.0 1134
NDIV 272 385
Hist update for nnode count 10464
d 1.0 198
d 3.0 56172
d 9.0 1134
NDIV 273 35
Hist update for nnode count 10496
d 1.0 194
d 3.0 56376
d 9.0 1134
NDIV 274 286
Hist update for nnode count 10516
d 1.0 194
d 3.0 56520
d 9.0 1134
NDIV 275 138
Hist update for nnode count 10532
d 1.0 194
d 3.0 56640
d 9.0 1134
NDIV 276 119
Hist update for nnode count 10548
d 1.0 194
d 3.0 56760
d 9.0 1134
NDIV 277 306
Hist update for nnode count 10570
d 1.0 192
d 3.0 56910
d 9.0 1134
NDIV 278 378
Hist update for nnode count 10594
d 1.0 188
d 3.0 57066
d 9.0 1134
NDIV 279 418
Hist update for nnode count 10612
d 1.0 186
d 3.0 57192
d 9.0 1134
NDIV 280 215
Hist update for nnode count 10624
d 1.0 186
d 3.0 57288
d 9.0 1134
NDIV 281 69
Hist update for nnode count 10656
d 1.0 178
d 3.0 57480
d 9.0 1134
NDIV 282 248
Hist update for nnode count 10702
d 1.0 164
d 3.0 57738
d 9.0 1134
NDIV 283 345
Hist update for nnode count 10740
d 1.0 154
d 3.0 57960
d 9.0 1134
NDIV 284 221
Hist update for nnode count 10774
d 1.0 148
d 3.0 58170
d 9.0 1134
NDIV 285 236
Hist update for nnode count 10794
d 1.0 148
d 3.0 58314
d 9.0 1134
NDIV 286 457
Hist update for nnode count 10816
d 1.0 146
d 3.0 58464
d 9.0 1134
NDIV 287 20
Hist update for nnode count 10838
d 1.0 144
d 3.0 58614
d 9.0 1134
NDIV 288 112
Hist update for nnode count 10854
d 1.0 144
d 3.0 58734
d 9.0 1134
NDIV 289 141
Hist update for nnode count 10872
d 1.0 142
d 3.0 58860
d 9.0 1134
NDIV 290 115
Hist update for nnode count 10880
d 1.0 142
d 3.0 58932
d 9.0 1134
NDIV 291 241
Hist update for nnode count 10892
d 1.0 142
d 3.0 59028
d 9.0 1134
NDIV 292 50
Hist update for nnode count 10918
d 1.0 140
d 3.0 59202
d 9.0 1134
NDIV 293 253
Hist update for nnode count 10936
d 1.0 138
d 3.0 59328
d 9.0 1134
NDIV 294 417
Hist update for nnode count 10962
d 1.0 136
d 3.0 59502
d 9.0 1134
NDIV 295 337
Hist update for nnode count 10992
d 1.0 130
d 3.0 59688
d 9.0 1134
NDIV 296 150
Hist update for nnode count 11014
d 1.0 128
d 3.0 59838
d 9.0 1134
NDIV 297 80
Hist update for nnode count 11034
d 1.0 128
d 3.0 59982
d 9.0 1134
NDIV 298 156
Hist update for nnode count 11058
d 1.0 124
d 3.0 60138
d 9.0 1134
NDIV 299 448
Hist update for nnode count 11100
d 1.0 114
d 3.0 60384
d 9.0 1134
Found alternative choice -- moving from  542 to 526
NDIV 300 542
Hist update for nnode count 11118
d 1.0 112
d 3.0 60510
d 9.0 1134
NDIV 301 214
Hist update for nnode count 11126
d 1.0 112
d 3.0 60582
d 9.0 1134
NDIV 302 376
Hist update for nnode count 11144
d 1.0 110
d 3.0 60708
d 9.0 1134
NDIV 303 311
Hist update for nnode count 11166
d 1.0 108
d 3.0 60858
d 9.0 1134
NDIV 304 66
Hist update for nnode count 11174
d 1.0 108
d 3.0 60930
d 9.0 1134
NDIV 305 393
Hist update for nnode count 11194
d 1.0 108
d 3.0 61074
d 9.0 1134
NDIV 306 308
Hist update for nnode count 11220
d 1.0 106
d 3.0 61248
d 9.0 1134
NDIV 307 85
Hist update for nnode count 11232
d 1.0 106
d 3.0 61344
d 9.0 1134
NDIV 308 468
Hist update for nnode count 11248
d 1.0 106
d 3.0 61464
d 9.0 1134
NDIV 309 321
Hist update for nnode count 11272
d 1.0 102
d 3.0 61620
d 9.0 1134
NDIV 310 243
Hist update for nnode count 11288
d 1.0 102
d 3.0 61740
d 9.0 1134
NDIV 311 373
Hist update for nnode count 11310
d 1.0 100
d 3.0 61890
d 9.0 1134
NDIV 312 235
Hist update for nnode count 11326
d 1.0 100
d 3.0 62010
d 9.0 1134
NDIV 313 329
Hist update for nnode count 11334
d 1.0 100
d 3.0 62082
d 9.0 1134
NDIV 314 277
Hist update for nnode count 11356
d 1.0 98
d 3.0 62232
d 9.0 1134
NDIV 315 469
Hist update for nnode count 11384
d 1.0 94
d 3.0 62412
d 9.0 1134
NDIV 316 380
Hist update for nnode count 11402
d 1.0 92
d 3.0 62538
d 9.0 1134
NDIV 317 268
Hist update for nnode count 11418
d 1.0 92
d 3.0 62658
d 9.0 1134
NDIV 318 53
Hist update for nnode count 11440
d 1.0 90
d 3.0 62808
d 9.0 1134
NDIV 319 220
Hist update for nnode count 11470
d 1.0 84
d 3.0 62994
d 9.0 1134
NDIV 320 262
Hist update for nnode count 11490
d 1.0 84
d 3.0 63138
d 9.0 1134
NDIV 321 232
Hist update for nnode count 11506
d 1.0 84
d 3.0 63258
d 9.0 1134
NDIV 322 510
Hist update for nnode count 11518
d 1.0 84
d 3.0 63354
d 9.0 1134
NDIV 323 227
Hist update for nnode count 11530
d 1.0 84
d 3.0 63450
d 9.0 1134
NDIV 324 252
Hist update for nnode count 11546
d 1.0 84
d 3.0 63570
d 9.0 1134
NDIV 325 48
Hist update for nnode count 11578
d 1.0 80
d 3.0 63774
d 9.0 1134
NDIV 326 108
Hist update for nnode count 11600
d 1.0 78
d 3.0 63924
d 9.0 1134
NDIV 327 288
Hist update for nnode count 11618
d 1.0 76
d 3.0 64050
d 9.0 1134
NDIV 328 123
Hist update for nnode count 11636
d 1.0 74
d 3.0 64176
d 9.0 1134
NDIV 329 105
Hist update for nnode count 11656
d 1.0 74
d 3.0 64320
d 9.0 1134
NDIV 330 82
Hist update for nnode count 11674
d 1.0 72
d 3.0 64446
d 9.0 1134
NDIV 331 504
Hist update for nnode count 11686
d 1.0 72
d 3.0 64542
d 9.0 1134
NDIV 332 283
Hist update for nnode count 11708
d 1.0 70
d 3.0 64692
d 9.0 1134
NDIV 333 328
Hist update for nnode count 11732
d 1.0 66
d 3.0 64848
d 9.0 1134
NDIV 334 314
Hist update for nnode count 11750
d 1.0 64
d 3.0 64974
d 9.0 1134
NDIV 335 399
Hist update for nnode count 11780
d 1.0 62
d 3.0 65172
d 9.0 1134
NDIV 336 165
Hist update for nnode count 11802
d 1.0 60
d 3.0 65322
d 9.0 1134
NDIV 337 223
Hist update for nnode count 11826
d 1.0 56
d 3.0 65478
d 9.0 1134
NDIV 338 144
Hist update for nnode count 11854
d 1.0 52
d 3.0 65658
d 9.0 1134
NDIV 339 503
Hist update for nnode count 11862
d 1.0 52
d 3.0 65730
d 9.0 1134
NDIV 340 307
Hist update for nnode count 11874
d 1.0 52
d 3.0 65826
d 9.0 1134
NDIV 341 87
Hist update for nnode count 11882
d 1.0 52
d 3.0 65898
d 9.0 1134
NDIV 342 91
Hist update for nnode count 11898
d 1.0 52
d 3.0 66018
d 9.0 1134
NDIV 343 449
Hist update for nnode count 11920
d 1.0 50
d 3.0 66168
d 9.0 1134
NDIV 344 374
Hist update for nnode count 11944
d 1.0 46
d 3.0 66324
d 9.0 1134
NDIV 345 139
Hist update for nnode count 11962
d 1.0 44
d 3.0 66450
d 9.0 1134
NDIV 346 287
Hist update for nnode count 11974
d 1.0 44
d 3.0 66546
d 9.0 1134
NDIV 347 459
Hist update for nnode count 11996
d 1.0 42
d 3.0 66696
d 9.0 1134
NDIV 348 104
Hist update for nnode count 12022
d 1.0 40
d 3.0 66870
d 9.0 1134
NDIV 349 213
Hist update for nnode count 12040
d 1.0 38
d 3.0 66996
d 9.0 1134
NDIV 350 348
Hist update for nnode count 12060
d 1.0 38
d 3.0 67140
d 9.0 1134
NDIV 351 339
Hist update for nnode count 12072
d 1.0 38
d 3.0 67236
d 9.0 1134
NDIV 352 330
Hist update for nnode count 12094
d 1.0 36
d 3.0 67386
d 9.0 1134
NDIV 353 192
Hist update for nnode count 12114
d 1.0 36
d 3.0 67530
d 9.0 1134
NDIV 354 352
Hist update for nnode count 12130
d 1.0 36
d 3.0 67650
d 9.0 1134
NDIV 355 428
Hist update for nnode count 12142
d 1.0 36
d 3.0 67746
d 9.0 1134
NDIV 356 383
Hist update for nnode count 12168
d 1.0 34
d 3.0 67920
d 9.0 1134
NDIV 357 292
Hist update for nnode count 12180
d 1.0 34
d 3.0 68016
d 9.0 1134
NDIV 358 147
Hist update for nnode count 12202
d 1.0 32
d 3.0 68166
d 9.0 1134
NDIV 359 369
Hist update for nnode count 12214
d 1.0 32
d 3.0 68262
d 9.0 1134
NDIV 360 498
Hist update for nnode count 12236
d 1.0 30
d 3.0 68412
d 9.0 1134
NDIV 361 359
Hist update for nnode count 12248
d 1.0 30
d 3.0 68508
d 9.0 1134
NDIV 362 377
Hist update for nnode count 12264
d 1.0 30
d 3.0 68628
d 9.0 1134
NDIV 363 486
Hist update for nnode count 12272
d 1.0 30
d 3.0 68700
d 9.0 1134
NDIV 364 261
Hist update for nnode count 12284
d 1.0 30
d 3.0 68796
d 9.0 1134
NDIV 365 56
Hist update for nnode count 12300
d 1.0 30
d 3.0 68916
d 9.0 1134
NDIV 366 482
Hist update for nnode count 12312
d 1.0 30
d 3.0 69012
d 9.0 1134
NDIV 367 174
Hist update for nnode count 12328
d 1.0 30
d 3.0 69132
d 9.0 1134
NDIV 368 382
Hist update for nnode count 12360
d 1.0 22
d 3.0 69324
d 9.0 1134
NDIV 369 260
Hist update for nnode count 12368
d 1.0 22
d 3.0 69396
d 9.0 1134
NDIV 370 325
Hist update for nnode count 12384
d 1.0 22
d 3.0 69516
d 9.0 1134
NDIV 371 143
Hist update for nnode count 12400
d 1.0 22
d 3.0 69636
d 9.0 1134
NDIV 372 534
Hist update for nnode count 12452
d 1.0 22
d 3.0 69616
d 9.0 1410
NDIV 373 487
Hist update for nnode count 12468
d 1.0 22
d 3.0 69736
d 9.0 1410
NDIV 374 327
Hist update for nnode count 12480
d 1.0 22
d 3.0 69832
d 9.0 1410
NDIV 375 415
Hist update for nnode count 12488
d 1.0 22
d 3.0 69904
d 9.0 1410
Found alternative choice -- moving from  541 to 11108
NDIV 376 541
Hist update for nnode count 12542
d 1.0 22
d 3.0 69882
d 9.0 1692
NDIV 377 490
Hist update for nnode count 12554
d 1.0 22
d 3.0 69978
d 9.0 1692
NDIV 378 478
Hist update for nnode count 12562
d 1.0 22
d 3.0 70050
d 9.0 1692
NDIV 379 22
Hist update for nnode count 12578
d 1.0 22
d 3.0 70170
d 9.0 1692
NDIV 380 296
Hist update for nnode count 12596
d 1.0 20
d 3.0 70296
d 9.0 1692
NDIV 381 313
Hist update for nnode count 12614
d 1.0 18
d 3.0 70422
d 9.0 1692
NDIV 382 175
Hist update for nnode count 12626
d 1.0 18
d 3.0 70518
d 9.0 1692
NDIV 383 137
Hist update for nnode count 12634
d 1.0 18
d 3.0 70590
d 9.0 1692
NDIV 384 489
Hist update for nnode count 12646
d 1.0 18
d 3.0 70686
d 9.0 1692
NDIV 385 297
Hist update for nnode count 12654
d 1.0 18
d 3.0 70758
d 9.0 1692
NDIV 386 361
Hist update for nnode count 12710
d 1.0 18
d 3.0 70734
d 9.0 1980
NDIV 387 100
Hist update for nnode count 12718
d 1.0 18
d 3.0 70806
d 9.0 1980
NDIV 388 433
Hist update for nnode count 12730
d 1.0 18
d 3.0 70902
d 9.0 1980
NDIV 389 532
Hist update for nnode count 12784
d 1.0 18
d 3.0 70880
d 9.0 2262
NDIV 390 78
Hist update for nnode count 12800
d 1.0 18
d 3.0 71000
d 9.0 2262
NDIV 391 316
Hist update for nnode count 12808
d 1.0 18
d 3.0 71072
d 9.0 2262
NDIV 392 120
Hist update for nnode count 12824
d 1.0 18
d 3.0 71192
d 9.0 2262
NDIV 393 544
Hist update for nnode count 12878
d 1.0 18
d 3.0 71170
d 9.0 2544
NDIV 394 15
Hist update for nnode count 12902
d 1.0 14
d 3.0 71326
d 9.0 2544
NDIV 395 73
Hist update for nnode count 12910
d 1.0 14
d 3.0 71398
d 9.0 2544
NDIV 396 4
Hist update for nnode count 12926
d 1.0 14
d 3.0 71518
d 9.0 2544
NDIV 397 294
Hist update for nnode count 12938
d 1.0 14
d 3.0 71614
d 9.0 2544
NDIV 398 237
Hist update for nnode count 12950
d 1.0 14
d 3.0 71710
d 9.0 2544
NDIV 399 499
Hist update for nnode count 12962
d 1.0 14
d 3.0 71806
d 9.0 2544
NDIV 400 445
Hist update for nnode count 12978
d 1.0 14
d 3.0 71926
d 9.0 2544
NDIV 401 128
Hist update for nnode count 12986
d 1.0 14
d 3.0 71998
d 9.0 2544
NDIV 402 279
Hist update for nnode count 13002
d 1.0 14
d 3.0 72118
d 9.0 2544
NDIV 403 177
Hist update for nnode count 13010
d 1.0 14
d 3.0 72190
d 9.0 2544
NDIV 404 229
Hist update for nnode count 13066
d 1.0 14
d 3.0 72166
d 9.0 2832
NDIV 405 443
Hist update for nnode count 13078
d 1.0 14
d 3.0 72262
d 9.0 2832
NDIV 406 473
Hist update for nnode count 13104
d 1.0 12
d 3.0 72436
d 9.0 2832
NDIV 407 363
Hist update for nnode count 13116
d 1.0 12
d 3.0 72532
d 9.0 2832
NDIV 408 188
Hist update for nnode count 13172
d 1.0 12
d 3.0 72508
d 9.0 3120
NDIV 409 46
Hist update for nnode count 13180
d 1.0 12
d 3.0 72580
d 9.0 3120
NDIV 410 13
Hist update for nnode count 13202
d 1.0 10
d 3.0 72730
d 9.0 3120
NDIV 411 77
Hist update for nnode count 13214
d 1.0 10
d 3.0 72826
d 9.0 3120
NDIV 412 525
Hist update for nnode count 13258
d 1.0 10
d 3.0 72810
d 9.0 3360
NDIV 413 30
Hist update for nnode count 13314
d 1.0 10
d 3.0 72786
d 9.0 3648
NDIV 414 496
Hist update for nnode count 13370
d 1.0 10
d 3.0 72762
d 9.0 3936
NDIV 415 96
Hist update for nnode count 13382
d 1.0 10
d 3.0 72858
d 9.0 3936
NDIV 416 463
Hist update for nnode count 13394
d 1.0 10
d 3.0 72954
d 9.0 3936
NDIV 417 233
Hist update for nnode count 13450
d 1.0 10
d 3.0 72930
d 9.0 4224
NDIV 418 300
Hist update for nnode count 13462
d 1.0 10
d 3.0 73026
d 9.0 4224
NDIV 419 509
Hist update for nnode count 13484
d 1.0 8
d 3.0 73176
d 9.0 4224
NDIV 420 205
Hist update for nnode count 13492
d 1.0 8
d 3.0 73248
d 9.0 4224
NDIV 421 162
Hist update for nnode count 13500
d 1.0 8
d 3.0 73320
d 9.0 4224
NDIV 422 461
Hist update for nnode count 13516
d 1.0 8
d 3.0 73440
d 9.0 4224
NDIV 423 160
Hist update for nnode count 13524
d 1.0 8
d 3.0 73512
d 9.0 4224
NDIV 424 391
Hist update for nnode count 13540
d 1.0 8
d 3.0 73632
d 9.0 4224
NDIV 425 196
Hist update for nnode count 13558
d 1.0 6
d 3.0 73758
d 9.0 4224
NDIV 426 88
Hist update for nnode count 13566
d 1.0 6
d 3.0 73830
d 9.0 4224
NDIV 427 255
Hist update for nnode count 13582
d 1.0 6
d 3.0 73950
d 9.0 4224
NDIV 428 146
Hist update for nnode count 13594
d 1.0 6
d 3.0 74046
d 9.0 4224
NDIV 429 212
Hist update for nnode count 13610
d 1.0 6
d 3.0 74166
d 9.0 4224
NDIV 430 170
Hist update for nnode count 13618
d 1.0 6
d 3.0 74238
d 9.0 4224
NDIV 431 26
Hist update for nnode count 13626
d 1.0 6
d 3.0 74310
d 9.0 4224
NDIV 432 333
Hist update for nnode count 13634
d 1.0 6
d 3.0 74382
d 9.0 4224
NDIV 433 400
Hist update for nnode count 13646
d 1.0 6
d 3.0 74478
d 9.0 4224
NDIV 434 3
Hist update for nnode count 13664
d 1.0 4
d 3.0 74604
d 9.0 4224
NDIV 435 183
Hist update for nnode count 13680
d 1.0 4
d 3.0 74724
d 9.0 4224
NDIV 436 149
Hist update for nnode count 13692
d 1.0 4
d 3.0 74820
d 9.0 4224
NDIV 437 155
Hist update for nnode count 13704
d 1.0 4
d 3.0 74916
d 9.0 4224
NDIV 438 485
Hist update for nnode count 13716
d 1.0 4
d 3.0 75012
d 9.0 4224
NDIV 439 264
Hist update for nnode count 13724
d 1.0 4
d 3.0 75084
d 9.0 4224
NDIV 440 338
Hist update for nnode count 13740
d 1.0 4
d 3.0 75204
d 9.0 4224
NDIV 441 209
Hist update for nnode count 13758
d 1.0 2
d 3.0 75330
d 9.0 4224
NDIV 442 522
Hist update for nnode count 13802
d 1.0 2
d 3.0 75314
d 9.0 4464
NDIV 443 173
Hist update for nnode count 13814
d 1.0 2
d 3.0 75410
d 9.0 4464
NDIV 444 419
Hist update for nnode count 13830
d 1.0 2
d 3.0 75530
d 9.0 4464
NDIV 445 386
Hist update for nnode count 13846
d 1.0 2
d 3.0 75650
d 9.0 4464
NDIV 446 79
Hist update for nnode count 13854
d 1.0 2
d 3.0 75722
d 9.0 4464
NDIV 447 65
Hist update for nnode count 13910
d 1.0 2
d 3.0 75698
d 9.0 4752
NDIV 448 202
Hist update for nnode count 13926
d 1.0 2
d 3.0 75818
d 9.0 4752
NDIV 449 500
Hist update for nnode count 13938
d 1.0 2
d 3.0 75914
d 9.0 4752
NDIV 450 452
Hist update for nnode count 13950
d 1.0 2
d 3.0 76010
d 9.0 4752
NDIV 451 83
Hist update for nnode count 13962
d 1.0 2
d 3.0 76106
d 9.0 4752
NDIV 452 387
Hist update for nnode count 13974
d 1.0 2
d 3.0 76202
d 9.0 4752
NDIV 453 31
Hist update for nnode count 14030
d 1.0 2
d 3.0 76178
d 9.0 5040
NDIV 454 375
Hist update for nnode count 14042
d 1.0 2
d 3.0 76274
d 9.0 5040
NDIV 455 154
Hist update for nnode count 14050
d 1.0 2
d 3.0 76346
d 9.0 5040
NDIV 456 273
Hist update for nnode count 14066
d 1.0 2
d 3.0 76466
d 9.0 5040
NDIV 457 515
Hist update for nnode count 14106
d 1.0 2
d 3.0 76454
d 9.0 5268
NDIV 458 379
Hist update for nnode count 14162
d 1.0 2
d 3.0 76430
d 9.0 5556
NDIV 459 204
Hist update for nnode count 14174
d 1.0 2
d 3.0 76526
d 9.0 5556
NDIV 460 203
Hist update for nnode count 14182
d 1.0 2
d 3.0 76598
d 9.0 5556
NDIV 461 136
Hist update for nnode count 14238
d 1.0 2
d 3.0 76574
d 9.0 5844
NDIV 462 7
Hist update for nnode count 14246
d 1.0 2
d 3.0 76646
d 9.0 5844
NDIV 463 481
Hist update for nnode count 14254
d 1.0 2
d 3.0 76718
d 9.0 5844
NDIV 464 389
Hist update for nnode count 14266
d 1.0 2
d 3.0 76814
d 9.0 5844
NDIV 465 421
Hist update for nnode count 14282
d 1.0 2
d 3.0 76934
d 9.0 5844
NDIV 466 257
Hist update for nnode count 14302
d 1.0 2
d 3.0 77078
d 9.0 5844
NDIV 467 436
Hist update for nnode count 14310
d 1.0 2
d 3.0 77150
d 9.0 5844
NDIV 468 530
Hist update for nnode count 14354
d 1.0 2
d 3.0 77134
d 9.0 6084
NDIV 469 71
Hist update for nnode count 14366
d 1.0 2
d 3.0 77230
d 9.0 6084
NDIV 470 275
Hist update for nnode count 14422
d 1.0 2
d 3.0 77206
d 9.0 6372
NDIV 471 513
Hist update for nnode count 14464
d 1.0 2
d 3.0 77192
d 9.0 6606
NDIV 472 395
Hist update for nnode count 14476
d 1.0 2
d 3.0 77288
d 9.0 6606
NDIV 473 61
Hist update for nnode count 14488
d 1.0 2
d 3.0 77384
d 9.0 6606
NDIV 474 113
Hist update for nnode count 14506
d 3.0 77510
d 9.0 6606
NDIV 475 270
Hist update for nnode count 14562
d 3.0 77486
d 9.0 6894
NDIV 476 551
Hist update for nnode count 14610
d 3.0 77470
d 9.0 7158
NDIV 477 432
Hist update for nnode count 14666
d 3.0 77446
d 9.0 7446
NDIV 478 508
Hist update for nnode count 14682
d 3.0 77566
d 9.0 7446
NDIV 479 12
Hist update for nnode count 14738
d 3.0 77542
d 9.0 7734
NDIV 480 346
Hist update for nnode count 14746
d 3.0 77614
d 9.0 7734
NDIV 481 403
Hist update for nnode count 14802
d 3.0 77590
d 9.0 8022
NDIV 482 408
Hist update for nnode count 14814
d 3.0 77686
d 9.0 8022
NDIV 483 52
Hist update for nnode count 14822
d 3.0 77758
d 9.0 8022
NDIV 484 29
Hist update for nnode count 14834
d 3.0 77854
d 9.0 8022
NDIV 485 210
Hist update for nnode count 14842
d 3.0 77926
d 9.0 8022
NDIV 486 494
Hist update for nnode count 14898
d 3.0 77902
d 9.0 8310
NDIV 487 444
Hist update for nnode count 14906
d 3.0 77974
d 9.0 8310
NDIV 488 505
Hist update for nnode count 14918
d 3.0 78070
d 9.0 8310
NDIV 489 42
Hist update for nnode count 14974
d 3.0 78046
d 9.0 8598
NDIV 490 312
Hist update for nnode count 14986
d 3.0 78142
d 9.0 8598
NDIV 491 446
Hist update for nnode count 14994
d 3.0 78214
d 9.0 8598
NDIV 492 2
Hist update for nnode count 15050
d 3.0 78190
d 9.0 8886
NDIV 493 319
Hist update for nnode count 15066
d 3.0 78310
d 9.0 8886
NDIV 494 410
Hist update for nnode count 15074
d 3.0 78382
d 9.0 8886
NDIV 495 548
Hist update for nnode count 15126
d 3.0 78362
d 9.0 9162
NDIV 496 109
Hist update for nnode count 15142
d 3.0 78482
d 9.0 9162
NDIV 497 282
Hist update for nnode count 15198
d 3.0 78458
d 9.0 9450
NDIV 498 8
Hist update for nnode count 15214
d 3.0 78578
d 9.0 9450
NDIV 499 126
Hist update for nnode count 15226
d 3.0 78674
d 9.0 9450
NDIV 500 55
Hist update for nnode count 15242
d 3.0 78794
d 9.0 9450
NDIV 501 304
Hist update for nnode count 15254
d 3.0 78890
d 9.0 9450
NDIV 502 365
Hist update for nnode count 15262
d 3.0 78962
d 9.0 9450
NDIV 503 553
Hist update for nnode count 15314
d 3.0 78942
d 9.0 9726
NDIV 504 39
Hist update for nnode count 15322
d 3.0 79014
d 9.0 9726
NDIV 505 540
Hist update for nnode count 15364
d 3.0 79000
d 9.0 9960
NDIV 506 514
Hist update for nnode count 15402
d 3.0 78990
d 9.0 10182
NDIV 507 480
Hist update for nnode count 15410
d 3.0 79062
d 9.0 10182
NDIV 508 244
Hist update for nnode count 15422
d 3.0 79158
d 9.0 10182
NDIV 509 506
Hist update for nnode count 15430
d 3.0 79230
d 9.0 10182
NDIV 510 49
Hist update for nnode count 15438
d 3.0 79302
d 9.0 10182
NDIV 511 281
Hist update for nnode count 15446
d 3.0 79374
d 9.0 10182
NDIV 512 23
Hist update for nnode count 15502
d 3.0 79350
d 9.0 10470
NDIV 513 559
Hist update for nnode count 15550
d 3.0 79334
d 9.0 10734
NDIV 514 81
Hist update for nnode count 15558
d 3.0 79406
d 9.0 10734
NDIV 515 554
Hist update for nnode count 15582
d 3.0 79402
d 9.0 10890
NDIV 516 435
Hist update for nnode count 15638
d 3.0 79378
d 9.0 11178
NDIV 517 484
Hist update for nnode count 15650
d 3.0 79474
d 9.0 11178
NDIV 518 536
Hist update for nnode count 15676
d 3.0 79468
d 9.0 11340
NDIV 519 416
Hist update for nnode count 15732
d 3.0 79444
d 9.0 11628
NDIV 520 265
Hist update for nnode count 15740
d 3.0 79516
d 9.0 11628
NDIV 521 118
Hist update for nnode count 15748
d 3.0 79588
d 9.0 11628
NDIV 522 62
Hist update for nnode count 15804
d 3.0 79564
d 9.0 11916
NDIV 523 464
Hist update for nnode count 15820
d 3.0 79684
d 9.0 11916
NDIV 524 467
Hist update for nnode count 15828
d 3.0 79756
d 9.0 11916
NDIV 525 201
Hist update for nnode count 15884
d 3.0 79732
d 9.0 12204
NDIV 526 169
Hist update for nnode count 15940
d 3.0 79708
d 9.0 12492
NDIV 527 455
Hist update for nnode count 15948
d 3.0 79780
d 9.0 12492
NDIV 528 564
Hist update for nnode count 15992
d 3.0 79764
d 9.0 12732
NDIV 529 228
Hist update for nnode count 16000
d 3.0 79836
d 9.0 12732
NDIV 530 462
Hist update for nnode count 16008
d 3.0 79908
d 9.0 12732
NDIV 531 247
Hist update for nnode count 16020
d 3.0 80004
d 9.0 12732
NDIV 532 394
Hist update for nnode count 16028
d 3.0 80076
d 9.0 12732
NDIV 533 412
Hist update for nnode count 16036
d 3.0 80148
d 9.0 12732
NDIV 534 41
Hist update for nnode count 16092
d 3.0 80124
d 9.0 13020
NDIV 535 251
Hist update for nnode count 16148
d 3.0 80100
d 9.0 13308
NDIV 536 299
Hist update for nnode count 16156
d 3.0 80172
d 9.0 13308
NDIV 537 272
Hist update for nnode count 16212
d 3.0 80148
d 9.0 13596
NDIV 538 305
Hist update for nnode count 16220
d 3.0 80220
d 9.0 13596
NDIV 539 211
Hist update for nnode count 16276
d 3.0 80196
d 9.0 13884
NDIV 540 358
Hist update for nnode count 16288
d 3.0 80292
d 9.0 13884
NDIV 541 7493
Hist update for nnode count 16344
d 3.0 80268
d 9.0 14172
NDIV 542 8049
Hist update for nnode count 16400
d 3.0 80244
d 9.0 14460
NDIV 543 5899
Hist update for nnode count 16444
d 3.0 80228
d 9.0 14700
NDIV 544 5713
Hist update for nnode count 16500
d 3.0 80204
d 9.0 14988
NDIV 545 8375
Hist update for nnode count 16544
d 3.0 80188
d 9.0 15228
NDIV 546 13587
Hist update for nnode count 16598
d 3.0 80166
d 9.0 15510
NDIV 547 14915
Hist update for nnode count 16654
d 3.0 80142
d 9.0 15798
NDIV 548 8562
Hist update for nnode count 16710
d 3.0 80118
d 9.0 16086
NDIV 549 13589
Hist update for nnode count 16764
d 3.0 80096
d 9.0 16368
NDIV 550 429
Hist update for nnode count 16820
d 3.0 80072
d 9.0 16656
NDIV 551 6843
Hist update for nnode count 16872
d 3.0 80052
d 9.0 16932
NDIV 552 7020
Hist update for nnode count 16928
d 3.0 80028
d 9.0 17220
NDIV 553 1538
Hist update for nnode count 16984
d 3.0 80004
d 9.0 17508
NDIV 554 1390
Hist update for nnode count 17040
d 3.0 79980
d 9.0 17796
NDIV 555 3583
Hist update for nnode count 17094
d 3.0 79958
d 9.0 18078
NDIV 556 9062
Hist update for nnode count 17150
d 3.0 79934
d 9.0 18366
NDIV 557 1772
Hist update for nnode count 17206
d 3.0 79910
d 9.0 18654
NDIV 558 2360
Hist update for nnode count 17262
d 3.0 79886
d 9.0 18942
NDIV 559 5868
Hist update for nnode count 17318
d 3.0 79862
d 9.0 19230
NDIV 560 4136
Hist update for nnode count 17374
d 3.0 79838
d 9.0 19518
NDIV 561 10522
Hist update for nnode count 17430
d 3.0 79814
d 9.0 19806
NDIV 562 8759
Hist update for nnode count 17484
d 3.0 79792
d 9.0 20088
NDIV 563 4190
Hist update for nnode count 17540
d 3.0 79768
d 9.0 20376
NDIV 564 12251
Hist update for nnode count 17596
d 3.0 79744
d 9.0 20664
NDIV 565 5774
Hist update for nnode count 17650
d 3.0 79722
d 9.0 20946
NDIV 566 5666
Hist update for nnode count 17706
d 3.0 79698
d 9.0 21234
NDIV 567 9112
Hist update for nnode count 17760
d 3.0 79676
d 9.0 21516
NDIV 568 142
Hist update for nnode count 17816
d 3.0 79652
d 9.0 21804
NDIV 569 11661
Hist update for nnode count 17870
d 3.0 79630
d 9.0 22086
NDIV 570 8017
Hist update for nnode count 17926
d 3.0 79606
d 9.0 22374
NDIV 571 13562
Hist update for nnode count 17982
d 3.0 79582
d 9.0 22662
NDIV 572 1241
Hist update for nnode count 18038
d 3.0 79558
d 9.0 22950
NDIV 573 3551
Hist update for nnode count 18092
d 3.0 79536
d 9.0 23232
NDIV 574 5245
Hist update for nnode count 18148
d 3.0 79512
d 9.0 23520
NDIV 575 9501
Hist update for nnode count 18204
d 3.0 79488
d 9.0 23808
NDIV 576 2823
Hist update for nnode count 18260
d 3.0 79464
d 9.0 24096
NDIV 577 1836
Hist update for nnode count 18304
d 3.0 79448
d 9.0 24336
NDIV 578 536
Hist update for nnode count 18348
d 3.0 79432
d 9.0 24576
NDIV 579 8540
Hist update for nnode count 18404
d 3.0 79408
d 9.0 24864
NDIV 580 12170
Hist update for nnode count 18460
d 3.0 79384
d 9.0 25152
NDIV 581 9454
Hist update for nnode count 18516
d 3.0 79360
d 9.0 25440
NDIV 582 5671
Hist update for nnode count 18570
d 3.0 79338
d 9.0 25722
NDIV 583 5640
Hist update for nnode count 18626
d 3.0 79314
d 9.0 26010
NDIV 584 4547
Hist update for nnode count 18682
d 3.0 79290
d 9.0 26298
NDIV 585 8614
Hist update for nnode count 18738
d 3.0 79266
d 9.0 26586
NDIV 586 627
Hist update for nnode count 18794
d 3.0 79242
d 9.0 26874
NDIV 587 7369
Hist update for nnode count 18838
d 3.0 79226
d 9.0 27114
NDIV 588 322
Hist update for nnode count 18846
d 3.0 79298
d 9.0 27114
NDIV 589 8283
Hist update for nnode count 18902
d 3.0 79274
d 9.0 27402
NDIV 590 2422
Hist update for nnode count 18956
d 3.0 79252
d 9.0 27684
NDIV 591 14269
Hist update for nnode count 19012
d 3.0 79228
d 9.0 27972
NDIV 592 8446
Hist update for nnode count 19066
d 3.0 79206
d 9.0 28254
NDIV 593 9115
Hist update for nnode count 19100
d 3.0 79196
d 9.0 28452
NDIV 594 3142
Hist update for nnode count 19156
d 3.0 79172
d 9.0 28740
NDIV 595 3093
Hist update for nnode count 19212
d 3.0 79148
d 9.0 29028
NDIV 596 4324
Hist update for nnode count 19268
d 3.0 79124
d 9.0 29316
NDIV 597 11572
Hist update for nnode count 19322
d 3.0 79102
d 9.0 29598
NDIV 598 5998
Hist update for nnode count 19376
d 3.0 79080
d 9.0 29880
NDIV 599 5164
Hist update for nnode count 19432
d 3.0 79056
d 9.0 30168
NDIV 600 3562
Hist update for nnode count 19486
d 3.0 79034
d 9.0 30450
NDIV 601 13931
Hist update for nnode count 19542
d 3.0 79010
d 9.0 30738
NDIV 602 2887
Hist update for nnode count 19598
d 3.0 78986
d 9.0 31026
NDIV 603 3640
Hist update for nnode count 19654
d 3.0 78962
d 9.0 31314
NDIV 604 2406
Hist update for nnode count 19710
d 3.0 78938
d 9.0 31602
NDIV 605 11080
Hist update for nnode count 19766
d 3.0 78914
d 9.0 31890
NDIV 606 8792
Hist update for nnode count 19822
d 3.0 78890
d 9.0 32178
NDIV 607 2142
Hist update for nnode count 19876
d 3.0 78868
d 9.0 32460
NDIV 608 4838
Hist update for nnode count 19932
d 3.0 78844
d 9.0 32748
NDIV 609 3884
Hist update for nnode count 19986
d 3.0 78822
d 9.0 33030
NDIV 610 1926
Hist update for nnode count 20042
d 3.0 78798
d 9.0 33318
NDIV 611 6282
Hist update for nnode count 20098
d 3.0 78774
d 9.0 33606
NDIV 612 7388
Hist update for nnode count 20154
d 3.0 78750
d 9.0 33894
NDIV 613 12038
Hist update for nnode count 20210
d 3.0 78726
d 9.0 34182
NDIV 614 3147
Hist update for nnode count 20264
d 3.0 78704
d 9.0 34464
NDIV 615 6351
Hist update for nnode count 20320
d 3.0 78680
d 9.0 34752
NDIV 616 2688
Hist update for nnode count 20376
d 3.0 78656
d 9.0 35040
NDIV 617 12142
Hist update for nnode count 20430
d 3.0 78634
d 9.0 35322
NDIV 618 10085
Hist update for nnode count 20486
d 3.0 78610
d 9.0 35610
NDIV 619 10564
Hist update for nnode count 20540
d 3.0 78588
d 9.0 35892
NDIV 620 6397
Hist update for nnode count 20596
d 3.0 78564
d 9.0 36180
NDIV 621 13279
Hist update for nnode count 20652
d 3.0 78564
d 9.0 36156
d 27.0 288
NDIV 622 13808
Hist update for nnode count 20708
d 3.0 78540
d 9.0 36444
d 27.0 288
NDIV 623 7934
Hist update for nnode count 20764
d 3.0 78516
d 9.0 36732
d 27.0 288
NDIV 624 683
Hist update for nnode count 20820
d 3.0 78492
d 9.0 37020
d 27.0 288
NDIV 625 3481
Hist update for nnode count 20876
d 3.0 78468
d 9.0 37308
d 27.0 288
NDIV 626 4107
Hist update for nnode count 20932
d 3.0 78444
d 9.0 37596
d 27.0 288
NDIV 627 9209
Hist update for nnode count 20988
d 3.0 78420
d 9.0 37884
d 27.0 288
NDIV 628 14977
Hist update for nnode count 21044
d 3.0 78396
d 9.0 38172
d 27.0 288
NDIV 629 1164
Hist update for nnode count 21100
d 3.0 78372
d 9.0 38460
d 27.0 288
NDIV 630 13619
Hist update for nnode count 21156
d 3.0 78348
d 9.0 38748
d 27.0 288
NDIV 631 12319
Hist update for nnode count 21212
d 3.0 78324
d 9.0 39036
d 27.0 288
NDIV 632 7664
Hist update for nnode count 21266
d 3.0 78302
d 9.0 39318
d 27.0 288
NDIV 633 11769
Hist update for nnode count 21322
d 3.0 78278
d 9.0 39606
d 27.0 288
NDIV 634 13928
Hist update for nnode count 21376
d 3.0 78256
d 9.0 39888
d 27.0 288
NDIV 635 3406
Hist update for nnode count 21432
d 3.0 78232
d 9.0 40176
d 27.0 288
NDIV 636 11154
Hist update for nnode count 21488
d 3.0 78208
d 9.0 40464
d 27.0 288
NDIV 637 6222
Hist update for nnode count 21544
d 3.0 78184
d 9.0 40752
d 27.0 288
NDIV 638 7289
Hist update for nnode count 21600
d 3.0 78160
d 9.0 41040
d 27.0 288
NDIV 639 1306
Hist update for nnode count 21656
d 3.0 78136
d 9.0 41328
d 27.0 288
NDIV 640 7712
Hist update for nnode count 21710
d 3.0 78114
d 9.0 41610
d 27.0 288
NDIV 641 11268
Hist update for nnode count 21766
d 3.0 78090
d 9.0 41898
d 27.0 288
NDIV 642 11682
Hist update for nnode count 21820
d 3.0 78068
d 9.0 42180
d 27.0 288
NDIV 643 4247
Hist update for nnode count 21864
d 3.0 78052
d 9.0 42420
d 27.0 288
NDIV 644 12214
Hist update for nnode count 21908
d 3.0 78036
d 9.0 42660
d 27.0 288
NDIV 645 1362
Hist update for nnode count 21964
d 3.0 78012
d 9.0 42948
d 27.0 288
NDIV 646 1170
Hist update for nnode count 22020
d 3.0 77988
d 9.0 43236
d 27.0 288
NDIV 647 10406
Hist update for nnode count 22076
d 3.0 77964
d 9.0 43524
d 27.0 288
NDIV 648 4585
Hist update for nnode count 22132
d 3.0 77940
d 9.0 43812
d 27.0 288
NDIV 649 7318
Hist update for nnode count 22188
d 3.0 77916
d 9.0 44100
d 27.0 288
NDIV 650 6173
Hist update for nnode count 22244
d 3.0 77892
d 9.0 44388
d 27.0 288
NDIV 651 8618
Hist update for nnode count 22288
d 3.0 77876
d 9.0 44628
d 27.0 288
NDIV 652 773
Hist update for nnode count 22344
d 3.0 77852
d 9.0 44916
d 27.0 288
NDIV 653 1800
Hist update for nnode count 22400
d 3.0 77828
d 9.0 45204
d 27.0 288
NDIV 654 13964
Hist update for nnode count 22456
d 3.0 77804
d 9.0 45492
d 27.0 288
NDIV 655 5594
Hist update for nnode count 22512
d 3.0 77780
d 9.0 45780
d 27.0 288
NDIV 656 2727
Hist update for nnode count 22566
d 3.0 77758
d 9.0 46062
d 27.0 288
NDIV 657 7418
Hist update for nnode count 22620
d 3.0 77736
d 9.0 46344
d 27.0 288
NDIV 658 12206
Hist update for nnode count 22674
d 3.0 77714
d 9.0 46626
d 27.0 288
NDIV 659 11506
Hist update for nnode count 22730
d 3.0 77690
d 9.0 46914
d 27.0 288
NDIV 660 10740
Hist update for nnode count 22784
d 3.0 77668
d 9.0 47196
d 27.0 288
NDIV 661 1907
Hist update for nnode count 22828
d 3.0 77652
d 9.0 47436
d 27.0 288
NDIV 662 11474
Hist update for nnode count 22884
d 3.0 77628
d 9.0 47724
d 27.0 288
NDIV 663 5434
Hist update for nnode count 22928
d 3.0 77612
d 9.0 47964
d 27.0 288
NDIV 664 6713
Hist update for nnode count 22982
d 3.0 77590
d 9.0 48246
d 27.0 288
NDIV 665 2555
Hist update for nnode count 23038
d 3.0 77566
d 9.0 48534
d 27.0 288
NDIV 666 8141
Hist update for nnode count 23082
d 3.0 77550
d 9.0 48774
d 27.0 288
NDIV 667 10617
Hist update for nnode count 23138
d 3.0 77526
d 9.0 49062
d 27.0 288
NDIV 668 7906
Hist update for nnode count 23192
d 3.0 77504
d 9.0 49344
d 27.0 288
NDIV 669 7300
Hist update for nnode count 23236
d 3.0 77488
d 9.0 49584
d 27.0 288
NDIV 670 2997
Hist update for nnode count 23292
d 3.0 77464
d 9.0 49872
d 27.0 288
NDIV 671 819
Hist update for nnode count 23346
d 3.0 77442
d 9.0 50154
d 27.0 288
NDIV 672 16018
Hist update for nnode count 23402
d 3.0 77418
d 9.0 50442
d 27.0 288
NDIV 673 8517
Hist update for nnode count 23458
d 3.0 77394
d 9.0 50730
d 27.0 288
NDIV 674 4837
Hist update for nnode count 23512
d 3.0 77372
d 9.0 51012
d 27.0 288
NDIV 675 7871
Hist update for nnode count 23564
d 3.0 77352
d 9.0 51288
d 27.0 288
NDIV 676 11330
Hist update for nnode count 23620
d 3.0 77328
d 9.0 51576
d 27.0 288
NDIV 677 11889
Hist update for nnode count 23676
d 3.0 77304
d 9.0 51864
d 27.0 288
NDIV 678 11083
Hist update for nnode count 23732
d 3.0 77280
d 9.0 52152
d 27.0 288
NDIV 679 3595
Hist update for nnode count 23784
d 3.0 77260
d 9.0 52428
d 27.0 288
NDIV 680 10421
Hist update for nnode count 23840
d 3.0 77236
d 9.0 52716
d 27.0 288
NDIV 681 10175
Hist update for nnode count 23896
d 3.0 77212
d 9.0 53004
d 27.0 288
NDIV 682 2677
Hist update for nnode count 23950
d 3.0 77190
d 9.0 53286
d 27.0 288
NDIV 683 12231
Hist update for nnode count 24006
d 3.0 77166
d 9.0 53574
d 27.0 288
NDIV 684 5556
Hist update for nnode count 24062
d 3.0 77142
d 9.0 53862
d 27.0 288
NDIV 685 7393
Hist update for nnode count 24116
d 3.0 77120
d 9.0 54144
d 27.0 288
NDIV 686 1723
Hist update for nnode count 24172
d 3.0 77096
d 9.0 54432
d 27.0 288
NDIV 687 2093
Hist update for nnode count 24228
d 3.0 77072
d 9.0 54720
d 27.0 288
NDIV 688 6616
Hist update for nnode count 24284
d 3.0 77048
d 9.0 55008
d 27.0 288
NDIV 689 6608
Hist update for nnode count 24336
d 3.0 77028
d 9.0 55284
d 27.0 288
NDIV 690 5064
Hist update for nnode count 24392
d 3.0 77004
d 9.0 55572
d 27.0 288
NDIV 691 13836
Hist update for nnode count 24448
d 3.0 76980
d 9.0 55860
d 27.0 288
NDIV 692 8589
Hist update for nnode count 24502
d 3.0 76958
d 9.0 56142
d 27.0 288
NDIV 693 10612
Hist update for nnode count 24558
d 3.0 76934
d 9.0 56430
d 27.0 288
NDIV 694 7599
Hist update for nnode count 24614
d 3.0 76910
d 9.0 56718
d 27.0 288
NDIV 695 9972
Hist update for nnode count 24670
d 3.0 76886
d 9.0 57006
d 27.0 288
NDIV 696 4918
Hist update for nnode count 24726
d 3.0 76862
d 9.0 57294
d 27.0 288
NDIV 697 10621
Hist update for nnode count 24780
d 3.0 76840
d 9.0 57576
d 27.0 288
NDIV 698 15258
Hist update for nnode count 24836
d 3.0 76816
d 9.0 57864
d 27.0 288
NDIV 699 8537
Hist update for nnode count 24892
d 3.0 76792
d 9.0 58152
d 27.0 288
NDIV 700 2932
Hist update for nnode count 24948
d 3.0 76768
d 9.0 58440
d 27.0 288
NDIV 701 3521
Hist update for nnode count 25002
d 3.0 76746
d 9.0 58722
d 27.0 288
NDIV 702 11093
Hist update for nnode count 25036
d 3.0 76736
d 9.0 58920
d 27.0 288
NDIV 703 5028
Hist update for nnode count 25092
d 3.0 76712
d 9.0 59208
d 27.0 288
NDIV 704 11221
Hist update for nnode count 25148
d 3.0 76688
d 9.0 59496
d 27.0 288
NDIV 705 15129
Hist update for nnode count 25192
d 3.0 76672
d 9.0 59736
d 27.0 288
NDIV 706 10061
Hist update for nnode count 25246
d 3.0 76650
d 9.0 60018
d 27.0 288
NDIV 707 6416
Hist update for nnode count 25300
d 3.0 76628
d 9.0 60300
d 27.0 288
NDIV 708 7227
Hist update for nnode count 25354
d 3.0 76606
d 9.0 60582
d 27.0 288
NDIV 709 10059
Hist update for nnode count 25408
d 3.0 76584
d 9.0 60864
d 27.0 288
NDIV 710 2131
Hist update for nnode count 25452
d 3.0 76568
d 9.0 61104
d 27.0 288
NDIV 711 8032
Hist update for nnode count 25506
d 3.0 76546
d 9.0 61386
d 27.0 288
NDIV 712 12814
Hist update for nnode count 25560
d 3.0 76524
d 9.0 61668
d 27.0 288
NDIV 713 11704
Hist update for nnode count 25616
d 3.0 76500
d 9.0 61956
d 27.0 288
NDIV 714 11463
Hist update for nnode count 25672
d 3.0 76476
d 9.0 62244
d 27.0 288
Found alternative choice -- moving from  12763 to 15577
NDIV 715 12763
Hist update for nnode count 25706
d 3.0 76466
d 9.0 62442
d 27.0 288
NDIV 716 2540
Hist update for nnode count 25762
d 3.0 76442
d 9.0 62730
d 27.0 288
NDIV 717 13626
Hist update for nnode count 25816
d 3.0 76420
d 9.0 63012
d 27.0 288
NDIV 718 5123
Hist update for nnode count 25870
d 3.0 76398
d 9.0 63294
d 27.0 288
NDIV 719 8955
Hist update for nnode count 25926
d 3.0 76374
d 9.0 63582
d 27.0 288
NDIV 720 6804
Hist update for nnode count 25982
d 3.0 76350
d 9.0 63870
d 27.0 288
NDIV 721 5295
Hist update for nnode count 26036
d 3.0 76328
d 9.0 64152
d 27.0 288
NDIV 722 13028
Hist update for nnode count 26092
d 3.0 76328
d 9.0 64128
d 27.0 576
NDIV 723 2331
Hist update for nnode count 26148
d 3.0 76304
d 9.0 64416
d 27.0 576
NDIV 724 13385
Hist update for nnode count 26204
d 3.0 76280
d 9.0 64704
d 27.0 576
NDIV 725 2908
Hist update for nnode count 26260
d 3.0 76256
d 9.0 64992
d 27.0 576
NDIV 726 4810
Hist update for nnode count 26316
d 3.0 76232
d 9.0 65280
d 27.0 576
NDIV 727 2148
Hist update for nnode count 26372
d 3.0 76208
d 9.0 65568
d 27.0 576
NDIV 728 7405
Hist update for nnode count 26428
d 3.0 76184
d 9.0 65856
d 27.0 576
NDIV 729 1784
Hist update for nnode count 26482
d 3.0 76162
d 9.0 66138
d 27.0 576
NDIV 730 14703
Hist update for nnode count 26538
d 3.0 76162
d 9.0 66114
d 27.0 864
NDIV 731 2819
Hist update for nnode count 26582
d 3.0 76146
d 9.0 66354
d 27.0 864
NDIV 732 7990
Hist update for nnode count 26636
d 3.0 76124
d 9.0 66636
d 27.0 864
Found alternative choice -- moving from  12428 to 15576
NDIV 733 12428
Hist update for nnode count 26678
d 3.0 76110
d 9.0 66870
d 27.0 864
NDIV 734 1666
Hist update for nnode count 26734
d 3.0 76086
d 9.0 67158
d 27.0 864
NDIV 735 3514
Hist update for nnode count 26788
d 3.0 76064
d 9.0 67440
d 27.0 864
NDIV 736 1235
Hist update for nnode count 26844
d 3.0 76040
d 9.0 67728
d 27.0 864
NDIV 737 817
Hist update for nnode count 26898
d 3.0 76018
d 9.0 68010
d 27.0 864
Found alternative choice -- moving from  13032 to 13017
NDIV 738 13032
Hist update for nnode count 26942
d 3.0 76002
d 9.0 68250
d 27.0 864
NDIV 739 5557
Hist update for nnode count 26994
d 3.0 75982
d 9.0 68526
d 27.0 864
NDIV 740 7067
Hist update for nnode count 27050
d 3.0 75958
d 9.0 68814
d 27.0 864
NDIV 741 11282
Hist update for nnode count 27106
d 3.0 75934
d 9.0 69102
d 27.0 864
NDIV 742 13079
Hist update for nnode count 27162
d 3.0 75910
d 9.0 69390
d 27.0 864
NDIV 743 1394
Hist update for nnode count 27216
d 3.0 75888
d 9.0 69672
d 27.0 864
NDIV 744 1029
Hist update for nnode count 27272
d 3.0 75864
d 9.0 69960
d 27.0 864
NDIV 745 11014
Hist update for nnode count 27328
d 3.0 75840
d 9.0 70248
d 27.0 864
NDIV 746 1167
Hist update for nnode count 27370
d 3.0 75826
d 9.0 70482
d 27.0 864
NDIV 747 11979
Hist update for nnode count 27426
d 3.0 75802
d 9.0 70770
d 27.0 864
NDIV 748 10233
Hist update for nnode count 27482
d 3.0 75778
d 9.0 71058
d 27.0 864
NDIV 749 4960
Hist update for nnode count 27526
d 3.0 75762
d 9.0 71298
d 27.0 864
NDIV 750 1968
Hist update for nnode count 27582
d 3.0 75738
d 9.0 71586
d 27.0 864
NDIV 751 11509
Hist update for nnode count 27636
d 3.0 75716
d 9.0 71868
d 27.0 864
NDIV 752 14086
Hist update for nnode count 27692
d 3.0 75716
d 9.0 71844
d 27.0 1152
NDIV 753 2879
Hist update for nnode count 27748
d 3.0 75692
d 9.0 72132
d 27.0 1152
NDIV 754 8648
Hist update for nnode count 27802
d 3.0 75670
d 9.0 72414
d 27.0 1152
NDIV 755 3545
Hist update for nnode count 27846
d 3.0 75654
d 9.0 72654
d 27.0 1152
NDIV 756 5632
Hist update for nnode count 27902
d 3.0 75630
d 9.0 72942
d 27.0 1152
NDIV 757 3605
Hist update for nnode count 27946
d 3.0 75614
d 9.0 73182
d 27.0 1152
NDIV 758 11053
Hist update for nnode count 27990
d 3.0 75598
d 9.0 73422
d 27.0 1152
NDIV 759 15231
Hist update for nnode count 28044
d 3.0 75576
d 9.0 73704
d 27.0 1152
NDIV 760 9570
Hist update for nnode count 28098
d 3.0 75554
d 9.0 73986
d 27.0 1152
NDIV 761 1491
Hist update for nnode count 28154
d 3.0 75530
d 9.0 74274
d 27.0 1152
NDIV 762 11849
Hist update for nnode count 28194
d 3.0 75518
d 9.0 74502
d 27.0 1152
NDIV 763 1787
Hist update for nnode count 28248
d 3.0 75496
d 9.0 74784
d 27.0 1152
NDIV 764 11946
Hist update for nnode count 28302
d 3.0 75474
d 9.0 75066
d 27.0 1152
"""