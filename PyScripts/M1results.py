import numpy as np
# Z = np.array([[[[-2.11356932e-11,  5.48450449e-10, -5.28076874e-09,
#                  2.23485971e-08, -8.80427443e-08],
#                 [ 1.04000119e-09, -2.69885806e-08,  2.59882697e-07,
#                   -1.09997544e-06,  4.20799743e-06],
#                 [-1.90555694e-08,  4.94525751e-07, -4.76231647e-06,
#                  2.01590301e-05, -7.52850260e-05],
#                 [ 1.54101473e-07, -3.99936285e-06,  3.85165596e-05,
#                   -1.63056765e-04,  5.97552433e-04],
#                 [-4.64142234e-07,  1.20461722e-05, -1.16019004e-04,
#                  4.91196672e-04, -1.77496900e-03]],
#
#                [[ 2.54378972e-08, -6.59966359e-07,  6.35365776e-06,
#                   -2.68874898e-05,  1.10128062e-04],
#                 [-1.25195285e-06,  3.24828292e-05, -3.12747161e-04,
#                  1.32364666e-03, -5.23982833e-03],
#                 [ 2.29429717e-05, -5.95300339e-04,  5.73201652e-03,
#                   -2.42622708e-02,  9.33099615e-02],
#                 [-1.85564634e-04,  4.81502335e-03, -4.63656833e-02,
#                  1.96272775e-01, -7.37175170e-01],
#                 [ 5.58971458e-04, -1.45046209e-02,  1.39677858e-01,
#                   -5.91323664e-01,  2.17983171e+00]],
#
#                [[-1.12353922e-05,  2.91426930e-04, -2.80512818e-03,
#                  1.18693809e-02, -5.18778320e-02],
#                 [ 5.53114748e-04, -1.43477214e-02,  1.38116156e-01,
#                   -5.84484510e-01,  2.45712479e+00],
#                 [-1.01386077e-02,  2.63006983e-01, -2.53198599e+00,
#                  1.07160748e+01, -4.35555740e+01],
#                 [ 8.20180878e-02, -2.12772635e+00,  2.04850583e+01,
#                   -8.67064420e+01,  3.42469981e+02],
#                 [-2.47102355e-01,  6.41057299e+00, -6.17222076e+01,
#                  2.61270304e+02, -1.00777224e+03]],
#
#                [[ 2.15533375e-03, -5.58933417e-02,  5.37904867e-01,
#                   -2.27577116e+00,  1.10996543e+01],
#                 [-1.06142927e-01,  2.75273616e+00, -2.64941712e+01,
#                  1.12105890e+02, -5.20164802e+02],
#                 [ 1.94617587e+00, -5.04752189e+01,  4.85845122e+02,
#                   -2.05600055e+03,  9.15439197e+03],
#                 [-1.57479232e+01,  4.08448793e+02, -3.93174857e+03,
#                  1.66399442e+04, -7.15767885e+04],
#                 [ 4.74553595e+01, -1.23087571e+03,  1.18491254e+04,
#                   -5.01519248e+04,  2.09587458e+05]],
#
#                [[-1.51345927e-01,  3.92437953e+00, -3.77650067e+01,
#                  1.59775708e+02, -1.46045605e+03],
#                 [ 7.45627291e+00, -1.93353422e+02,  1.86085574e+03,
#                   -7.87391593e+03,  6.15534965e+04],
#                 [-1.36761391e+02,  3.54664392e+03, -3.41360923e+04,
#                  1.44457671e+05, -1.00244730e+06],
#                 [ 1.10696872e+03, -2.87083813e+04,  2.76334138e+05,
#                   -1.16950796e+06,  7.40850579e+06],
#                 [-3.33665272e+03,  8.65365901e+04, -8.33010426e+05,
#                  3.52577881e+06, -2.08139152e+07]]],
#
#
#               [[[ 7.93256307e-11, -2.05824749e-09,  1.98160143e-08,
#                   -8.38540305e-08,  3.28636141e-07],
#                 [-3.90314714e-09,  1.01280226e-07, -9.75172270e-07,
#                  4.12707065e-06, -1.56997124e-05],
#                 [ 7.15135512e-08, -1.85574885e-06,  1.78693299e-05,
#                   -7.56336337e-05,  2.80670434e-04],
#                 [-5.78308637e-07,  1.50074838e-05, -1.44518892e-04,
#                  6.11746682e-04, -2.22557719e-03],
#                 [ 1.74177275e-06, -4.52016143e-05,  4.35305719e-04,
#                   -1.84279335e-03,  6.60340594e-03]],
#
#                [[-9.54603678e-08,  2.47642940e-06, -2.38388838e-05,
#                  1.00870869e-04, -4.12067626e-04],
#                 [ 4.69802508e-06, -1.21883201e-04,  1.17338828e-03,
#                   -4.96562924e-03,  1.96015063e-02],
#                 [-8.60921781e-05,  2.23363990e-03, -2.15051865e-02,
#                  9.10167838e-02, -3.48869476e-01],
#                 [ 6.96300932e-04, -1.80660619e-02,  1.73948569e-01,
#                   -7.36273019e-01,  2.75395625e+00],
#                 [-2.09739439e-03,  5.44202431e-02, -5.24011816e-01,
#                  2.21816441e+00, -8.13530004e+00]],
#
#                [[ 4.21567679e-05, -1.09337675e-03,  1.05232405e-02,
#                   -4.45222210e-02,  1.94566488e-01],
#                 [-2.07530225e-03,  5.38282851e-02, -5.18118589e-01,
#                  2.19235133e+00, -9.21636249e+00],
#                 [ 3.80393076e-02, -9.86696140e-01,  9.49806070e+00,
#                   -4.01940964e+01,  1.63333850e+02],
#                 [-3.07718116e-01,  7.98217800e+00, -7.68423615e+01,
#                  3.25213190e+02, -1.28361493e+03],
#                 [ 9.27065288e-01, -2.40487666e+01,  2.31523783e+02,
#                   -9.79936322e+02,  3.77445660e+03]],
#
#                [[-8.08580129e-03,  2.09666141e-01, -2.01757088e+00,
#                  8.53497447e+00, -4.17264795e+01],
#                 [ 3.98188226e-01, -1.03257612e+01,  9.93718584e+01,
#                   -4.20428689e+02,  1.95578903e+03],
#                 [-7.30078002e+00,  1.89332830e+02, -1.82222327e+03,
#                  7.71042439e+03, -3.44168723e+04],
#                 [ 5.90746503e+01, -1.53206216e+03,  1.47462277e+04,
#                   -6.24020425e+04,  2.69019940e+05],
#                 [-1.78014175e+02,  4.61683921e+03, -4.44399574e+04,
#                  1.88073269e+05, -7.87352343e+05]],
#
#                [[ 5.67676176e-01, -1.47183535e+01,  1.41622211e+02,
#                   -5.99102650e+02,  5.59134873e+03],
#                 [-2.79667839e+01,  7.25155863e+02, -6.97824657e+03,
#                  2.95238941e+04, -2.35642088e+05],
#                 [ 5.12951036e+02, -1.33011499e+04,  1.28008847e+05,
#                   -5.41647281e+05,  3.83333410e+06],
#                 [-4.15183563e+03,  1.07664641e+05, -1.03622470e+06,
#                  4.38503256e+06, -2.82857007e+07],
#                 [ 1.25143655e+04, -3.24532014e+05,  3.12365888e+06,
#                   -1.32196174e+07,  7.93352469e+07]]],
#
#
#               [[[-1.11639436e-10,  2.89644093e-09, -2.78831693e-08,
#                  1.17978815e-07, -4.60280936e-07],
#                 [ 5.49291617e-09, -1.42520083e-07,  1.37211937e-06,
#                   -5.80640245e-06,  2.19802786e-05],
#                 [-1.00637914e-07,  2.61129647e-06, -2.51422860e-05,
#                  1.06406100e-04, -3.92689048e-04],
#                 [ 8.13803121e-07, -2.11169686e-05,  2.03333137e-04,
#                   -8.60618299e-04,  3.11107260e-03],
#                 [-2.45097344e-06,  6.36012078e-05, -6.12443384e-04,
#                  2.59241153e-03, -9.22104108e-03]],
#
#                [[ 1.34329560e-07, -3.48446821e-06,  3.35393150e-05,
#                   -1.41901558e-04,  5.78470448e-04],
#                 [-6.61072631e-06,  1.71490536e-04, -1.65080657e-03,
#                  6.98525578e-03, -2.75135067e-02],
#                 [ 1.21139037e-04, -3.14265276e-03,  3.02541515e-02,
#                   -1.28031619e-01,  4.89468713e-01],
#                 [-9.79726938e-04,  2.54176100e-02, -2.44709573e-01,
#                  1.03567465e+00, -3.86111380e+00],
#                 [ 2.95105156e-03, -7.65632957e-02,  7.37157726e-01,
#                   -3.12009468e+00,  1.13955131e+01]],
#
#                [[-5.93133404e-05,  1.53821054e-03, -1.48030712e-02,
#                  6.26225616e-02, -2.73728869e-01],
#                 [ 2.91980180e-03, -7.57257965e-02,  7.28818381e-01,
#                   -3.08356038e+00,  1.29687614e+01],
#                 [-5.35171141e-02,  1.38805023e+00, -1.33602303e+01,
#                  5.65319382e+01, -2.29804247e+02],
#                 [ 4.32914568e-01, -1.12287796e+01,  1.08086009e+02,
#                   -4.57393439e+02,  1.80526164e+03],
#                 [-1.30421559e+00,  3.38293927e+01, -3.25652917e+02,
#                  1.37819472e+03, -5.30497479e+03]],
#
#                [[ 1.13746194e-02, -2.94918206e-01,  2.83763748e+00,
#                   -1.20027379e+01,  5.88310447e+01],
#                 [-5.60132715e-01,  1.45239423e+01, -1.39759446e+02,
#                  5.91235664e+02, -2.75820653e+03],
#                 [ 1.02697889e+01, -2.66304496e+02,  2.56277175e+03,
#                   -1.08427045e+04,  4.85368341e+04],
#                 [-8.30967450e+01,  2.15486431e+03, -2.07386578e+04,
#                  8.77505709e+04, -3.79307337e+05],
#                 [ 2.50396817e+02, -6.49351712e+03,  6.24978979e+04,
#                   -2.64466505e+05,  1.10969971e+06]],
#
#                [[-7.98426578e-01,  2.06991065e+01, -1.99148495e+02,
#                  8.42355146e+02, -8.03013801e+03],
#                 [ 3.93339482e+01, -1.01979977e+03,  9.81259053e+03,
#                   -4.15106870e+04,  3.38444656e+05],
#                 [-7.21427183e+02,  1.87053029e+04, -1.79998895e+05,
#                  7.61545509e+05, -5.50031428e+06],
#                 [ 5.83914223e+03, -1.51405397e+05,  1.45705831e+06,
#                   -6.16517950e+06,  4.05273007e+07],
#                 [-1.75999195e+04,  4.56372122e+05, -4.39218157e+06,
#                  1.85859988e+07, -1.13489905e+08]]],
#
#
#               [[[ 6.98256565e-11, -1.81144382e-09,  1.74365618e-08,
#                   -7.37694835e-08,  2.86672115e-07],
#                 [-3.43545458e-09,  8.91292975e-08, -8.58015712e-07,
#                  3.63048713e-06, -1.36857474e-05],
#                 [ 6.29401843e-08, -1.63299964e-06,  1.57214996e-05,
#                   -6.65289306e-05,  2.44361989e-04],
#                 [-5.08946397e-07,  1.32052950e-05, -1.27140623e-04,
#                  5.38073877e-04, -1.93440652e-03],
#                 [ 1.53277591e-06, -3.97712789e-05,  3.82939431e-04,
#                   -1.62077811e-03,  5.72792043e-03]],
#
#                [[-8.40065026e-08,  2.17891201e-06, -2.09707931e-05,
#                  8.87157398e-05, -3.61078880e-04],
#                 [ 4.13405018e-06, -1.07233132e-04,  1.03215027e-03,
#                   -4.36699162e-03,  1.71732159e-02],
#                 [-7.57525335e-05,  1.96504307e-03, -1.89155491e-02,
#                  8.00396033e-02, -3.05406562e-01],
#                 [ 6.12640445e-04, -1.58927220e-02,  1.52993573e-01,
#                   -6.47440201e-01,  2.40769633e+00],
#                 [-1.84529566e-03,  4.78710600e-02, -4.60862979e-01,
#                  1.95044453e+00, -7.10016696e+00]],
#
#                [[ 3.70876821e-05, -9.61731317e-04,  9.25435031e-03,
#                   -3.91450007e-02,  1.71200252e-01],
#                 [-1.82564998e-03,  4.73444504e-02, -4.55618275e-01,
#                  1.92746339e+00, -8.11348456e+00],
#                 [ 3.34614536e-02, -8.67798239e-01,  8.35188551e+00,
#                   -3.53359500e+01,  1.43764379e+02],
#                 [-2.70671970e-01,  7.01997324e+00, -6.75662415e+01,
#                  2.85892607e+02, -1.12901492e+03],
#                 [ 8.15418048e-01, -2.11488749e+01,  2.03566245e+02,
#                   -8.61419010e+02,  3.31596108e+03]],
#
#                [[-7.11118842e-03,  1.84359762e-01, -1.77368308e+00,
#                  7.50151694e+00, -3.68689876e+01],
#                 [ 3.50174777e-01, -9.07900233e+00,  8.73554161e+01,
#                   -3.69504291e+02,  1.72910243e+03],
#                 [-6.42014845e+00,  1.66464700e+02, -1.60180293e+03,
#                  6.77621944e+03, -3.04293909e+04],
#                 [ 5.19467139e+01, -1.34695928e+03,  1.29619724e+04,
#                   -5.48392517e+04,  2.37768593e+05],
#                 [-1.56528765e+02,  4.05887945e+03, -3.90614041e+04,
#                  1.65274037e+05, -6.95400165e+05]],
#
#                [[ 4.99068213e-01, -1.29370225e+01,  1.24455115e+02,
#                   -5.26355220e+02,  5.12721585e+03],
#                 [-2.45857276e+01,  6.37365755e+02, -6.13212524e+03,
#                  2.59379608e+04, -2.16136673e+05],
#                 [ 4.50920011e+02, -1.16904299e+04,  1.12483698e+05,
#                   -4.75844116e+05,  3.50965426e+06],
#                 [-3.64962536e+03,  9.46236537e+04, -9.10520519e+05,
#                  3.85219359e+06, -2.58252269e+07],
#                 [ 1.10002557e+04, -2.85213902e+05,  2.74464838e+06,
#                   -1.16129441e+07,  7.22104697e+07]]],
#
#
#               [[[-1.63764429e-11,  4.24807364e-10, -4.08871016e-09,
#                  1.72964089e-08, -6.69882468e-08],
#                 [ 8.05698413e-10, -2.09012268e-08,  2.01189469e-07,
#                   -8.51194998e-07,  3.19737131e-06],
#                 [-1.47605100e-08,  3.82932995e-07, -3.68629214e-06,
#                  1.55977110e-05, -5.70619415e-05],
#                 [ 1.19352529e-07, -3.09650192e-06,  2.98103451e-05,
#                   -1.26147742e-04,  4.51389127e-04],
#                 [-3.59439287e-07,  9.32567748e-06, -8.97843529e-05,
#                  3.79970043e-04, -1.33541256e-03]],
#
#                [[ 1.96997585e-08, -5.10915843e-07,  4.91679131e-06,
#                   -2.07979325e-05,  8.45521403e-05],
#                 [-9.69413126e-07,  2.51434124e-05, -2.41989116e-04,
#                  1.02373672e-03, -4.02159901e-03],
#                 [ 1.77630158e-05, -4.60738021e-04,  4.43464642e-03,
#                   -1.87628330e-02,  7.15014276e-02],
#                 [-1.43652335e-04,  3.72621649e-03, -3.58675249e-02,
#                  1.51768520e-01, -5.63398179e-01],
#                 [ 4.32674657e-04, -1.12235855e-02,  1.08041138e-01,
#                   -4.57198806e-01,  1.66023182e+00]],
#
#                [[-8.69587116e-06,  2.25474596e-04, -2.16943014e-03,
#                  9.17544435e-03, -4.01616923e-02],
#                 [ 4.28043254e-04, -1.10994151e-02,  1.06804238e-01,
#                   -4.51777970e-01,  1.90404527e+00],
#                 [-7.84517834e-03,  2.03440794e-01, -1.95776503e+00,
#                  8.28218291e+00, -3.37400285e+01],
#                 [ 6.34585553e-02, -1.64567460e+00,  1.58378202e+01,
#                   -6.70071259e+01,  2.64912023e+02],
#                 [-1.91168790e-01,  4.95776425e+00, -4.77157611e+01,
#                  2.01893961e+02, -7.77714543e+02]],
#
#                [[ 1.66706625e-03, -4.32151301e-02,  4.15719109e-01,
#                   -1.75801423e+00,  8.66504947e+00],
#                 [-8.20888260e-02,  2.12812325e+00, -2.04740266e+01,
#                  8.65929995e+01, -4.06534421e+02],
#                 [ 1.50499062e+00, -3.90185165e+01,  3.75416118e+02,
#                   -1.58796766e+03,  7.15531111e+03],
#                 [-1.21769137e+01,  3.15714021e+02, -3.03784869e+03,
#                  1.28510132e+04, -5.59067399e+04],
#                 [ 3.66914197e+01, -9.51342951e+02,  9.15450398e+03,
#                   -3.87295817e+04,  1.63472780e+05]],
#
#                [[-1.16973989e-01,  3.03194475e+00, -2.91643556e+01,
#                  1.23329278e+02, -1.22798727e+03],
#                 [ 5.76239412e+00, -1.49371129e+02,  1.43695174e+03,
#                   -6.07736347e+03,  5.17814444e+04],
#                 [-1.05684377e+02,  2.73968260e+03, -2.63580377e+04,
#                  1.11490264e+05, -8.40242578e+05],
#                 [ 8.55365465e+02, -2.21749183e+04,  2.13356664e+05,
#                   -9.02555280e+05,  6.17524962e+06],
#                 [-2.57809546e+03,  6.68384172e+04, -6.43127001e+05,
#                  2.72083390e+06, -1.72422837e+07]]]])


Z = np.array([[[-5.29786276e-08,  2.48190138e-06, -4.36467396e-05,
                3.41613412e-04, -1.00388563e-03],
               [ 6.79416305e-05, -3.16270001e-03,  5.52311589e-02,
                 -4.29093841e-01,  1.25155414e+00],
               [-3.32555904e-02,  1.53996226e+00, -2.67377436e+01,
                2.06375916e+02, -5.97638900e+02],
               [ 7.52931603e+00, -3.44258958e+02,  5.92785082e+03,
                 -4.54599685e+04,  1.30864061e+05],
               [-1.20977736e+03,  4.91977143e+04, -7.75730597e+05,
                5.57280678e+06, -1.52791169e+07]],

              [[ 1.97087150e-07, -9.22418375e-06,  1.61981558e-04,
                 -1.26546635e-03,  3.71089550e-03],
               [-2.53819695e-04,  1.18101019e-02, -2.06037999e-01,
                1.59838915e+00, -4.65354967e+00],
               [ 1.24722518e-01, -5.77657075e+00,  1.00260536e+02,
                 -7.73220978e+02,  2.23636292e+03],
               [-2.83380561e+01,  1.29617306e+03, -2.23181437e+04,
                1.71090077e+05, -4.92169167e+05],
               [ 4.65151152e+03, -1.89318900e+05,  2.98336068e+06,
                 -2.14036564e+07,  5.85855190e+07]],

              [[-2.75218513e-07,  1.28708420e-05, -2.25728825e-04,
                1.76051812e-03, -5.15235364e-03],
               [ 3.55878800e-04, -1.65544380e-02,  2.88573325e-01,
                 -2.23582595e+00,  6.49859703e+00],
               [-1.75501890e-01,  8.13124093e+00, -1.41103748e+02,
                1.08750617e+03, -3.14203172e+03],
               [ 4.00052379e+01, -1.83072332e+03,  3.15251941e+04,
                 -2.41613709e+05,  6.94666202e+05],
               [-6.70886990e+03,  2.73322367e+05, -4.30541825e+06,
                3.08526268e+07, -8.43206170e+07]],

              [[ 1.70970002e-07, -7.99066464e-06,  1.39984183e-04,
                 -1.09011015e-03,  3.18445556e-03],
               [-2.21932639e-04,  1.03227004e-02, -1.79830277e-01,
                1.39178232e+00, -4.03933721e+00],
               [ 1.09806557e-01, -5.09002259e+00,  8.83276717e+01,
                 -6.80436495e+02,  1.96420539e+03],
               [-2.51046346e+01,  1.14952541e+03, -1.97991584e+04,
                1.51728156e+05, -4.36063030e+05],
               [ 4.30171302e+03, -1.75450196e+05,  2.76312980e+06,
                 -1.97809509e+07,  5.39871172e+07]],

              [[-3.98631988e-08,  1.86226868e-06, -3.25933297e-05,
                2.53471864e-04, -7.39197370e-04],
               [ 5.19354616e-05, -2.41584724e-03,  4.20673274e-02,
                 -3.25281939e-01,  9.42832391e-01],
               [-2.57729866e-02,  1.19546168e+00, -2.07481008e+01,
                1.59787272e+02, -4.60935831e+02],
               [ 5.90836772e+00, -2.70727997e+02,  4.66448549e+03,
                 -3.57465195e+04,  1.02708254e+05],
               [-1.03459045e+03,  4.22496735e+04, -6.65354195e+05,
                4.75927596e+06, -1.29732106e+07]]])