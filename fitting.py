import numpy as np
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt

#モデルの作成

#最初の日の出からの時間 t0
#最後の日の出からの時間 t1

#t0 から t1を計算
def t01(t0):
    t1 = t0%1
    return(t1)
#昼夜を決める関数
def Light(t0,dusk):
    t1 = t01(t0)
    if t1 < dusk:
        return(1)
    else:
        return(0)

#ベータL（昼用の分解速度）
def BetaL(F,a,r,k,C0,tL):
    F1 = t01(F)
    tD = 1 - tL
    if F1 >= tL:
        F1 = (tL/tD)*(1-F1)
    bbunbo = a*tD*F1 + C0
    if bbunbo <= 0:
        b = 100000
    else:
        b = a*(r-tD)/(bbunbo**k)
    return(b)
#ベータD（夜用の分解速度）
def BetaD(F,a,r,k,C0,tL):
    F1 = t01(F)
    tD = 1 - tL
    if F1 <= tL:
        F1 = 1-(tD/tL)*F1
    bbunbo = a*tL*(1-F1) + C0
    if bbunbo <= 0:
        b = 100000
    else:
        b = a*tL/(bbunbo**k)
    return(b)
#ベータまとめ
def Beta(F,t0,a,r,k,C0,tL,dusk):
    Lt = Light(t0,dusk)
    bL = BetaL(F,a,r,k,C0,tL)
    bD = BetaD(F,a,r,k,C0,tL)
    b = Lt*bL + (1-Lt)*bD
    return(b)
#ベータ動作確認
def test_Beta(a,r,k,C0,tL,dusk):
    Fs = np.linspace(0,3,300)
    ts = np.linspace(0,3,300)
    bs = [Beta(Fs[i],ts[i],a,r,k,C0,tL,dusk) for i in range(len(Fs))]
    plt.plot(Fs,bs)
    plt.show()

#Cについての微分方程式
def dCdt(C,F,t0,a,r,k,C0,tL,dusk,lamda):
    beta = Beta(F,t0,a,r,k,C0,tL,dusk)
    Lt = Light(t0,dusk)
    if C < 0:
        C = 0
    dcdt = a*r*Lt - lamda*beta*C**k
    return(dcdt)
#Sについての微分方程式
def dSdt(C,S,F,t0,a,r,k,C0,tL,dusk,lamda,H):
    beta = Beta(F,t0,a,r,k,C0,tL,dusk)
    Lt = Light(t0,dusk)
    if C < 0:
        C = 0
    dsdt = a*(1-r)*Lt + lamda*beta*C**k - H*S
    return(dsdt)
#Fについての微分方程式
def dFdt(C,S,F,t0,a,r,k,C0,tL,dusk,lamda,H,b,w,K,n):
    dS = dSdt(C,S,F,t0,a,r,k,C0,tL,dusk,lamda,H)
    sgn = np.sign(dS)
    absdS = np.abs(dS)
    hil = (absdS**n)/((K**n)+(absdS**n))
    fs = sgn*hil
    Zs = tL - b*F
    dfdt = w + Zs*fs
    return(dfdt)
#連立微分方程式モデル
def model(p,photoperiod):
    dt = 0.01
    T = 1
    N = int(T/dt)
    ts,Cs,Ss,Fs,Bs = [],[],[],[],[]
    #parameters
    a = p[0]
    r = p[1]
    k = p[2]
    C0_LD = p[3]
    C0_SD = p[4]
    tL = p[5]
    dusk_LD = p[6]
    dusk_SD = p[7]
    lamda = p[8]
    H = p[9]
    b = p[10]
    w = p[11]
    K = p[12]
    n = p[13]
    S0_LD = p[14]
    S0_SD = p[15]
    F0 = 0
    if photoperiod == "LD":
        dusk = dusk_LD
        C0 = C0_LD
        S0 = S0_LD
    else:
        dusk = dusk_SD
        C0 = C0_SD
        S0 = S0_SD
    C,S,F = C0,S0,F0
    for i in range(N):
        t0 = i*dt
        B = Beta(F,t0,a,r,k,C0,tL,dusk)
        ts.append(t0)
        Cs.append(C)
        Ss.append(S)
        Fs.append(F)
        Bs.append(B)
        kc1 = dt*dCdt(C,F,t0,a,r,k,C0,tL,dusk,lamda)
        ks1 = dt*dSdt(C,S,F,t0,a,r,k,C0,tL,dusk,lamda,H)
        kf1 = dt*dFdt(C,S,F,t0,a,r,k,C0,tL,dusk,lamda,H,b,w,K,n)
        kc2 = dt*dCdt(C+0.5*kc1,F+0.5*kf1,t0+0.5*dt,a,r,k,C0,tL,dusk,lamda)
        ks2 = dt*dSdt(C+0.5*kc1,S+0.5*ks1,F+0.5*kf1,t0+0.5*dt,a,r,k,C0,tL,dusk,lamda,H)
        kf2 = dt*dFdt(C+0.5*kc1,S+0.5*ks1,F+0.5*kf1,t0+0.5*dt,a,r,k,C0,tL,dusk,lamda,H,b,w,K,n)
        kc3 = dt*dCdt(C+0.5*kc2,F+0.5*kf2,t0+0.5*dt,a,r,k,C0,tL,dusk,lamda)
        ks3 = dt*dSdt(C+0.5*kc2,S+0.5*ks2,F+0.5*kf2,t0+0.5*dt,a,r,k,C0,tL,dusk,lamda,H)
        kf3 = dt*dFdt(C+0.5*kc2,S+0.5*ks2,F+0.5*kf2,t0+0.5*dt,a,r,k,C0,tL,dusk,lamda,H,b,w,K,n)
        kc4 = dt*dCdt(C+kc3,F+kf3,t0+dt,a,r,k,C0,tL,dusk,lamda)
        ks4 = dt*dSdt(C+kc3,S+ks3,F+kf3,t0+dt,a,r,k,C0,tL,dusk,lamda,H)
        kf4 = dt*dFdt(C+kc3,S+ks3,F+kf3,t0+dt,a,r,k,C0,tL,dusk,lamda,H,b,w,K,n)
        kc = (kc1+2*kc2+2*kc3+kc4)/6
        ks = (ks1+2*ks2+2*ks3+ks4)/6
        kf = (kf1+2*kf2+2*kf3+kf4)/6
        C += kc
        S += ks
        F += kf
        #rest
        #夜明けと日没にリセットする
        t1 = t01(t0)
        t1next = t01(t0+dt)
        if (t1 > t1next)  or ((t1-dusk)*(t1next-dusk)<0):
            F = t0
    M = pd.DataFrame()
    M["t"] = ts
    M["C"] = Cs
    M["S"] = Ss
    M["F"] = Fs
    M["B"] = Bs
    return(M)
#モデルの確認
def test_model():
    a = 144
    r = 0.6
    k = 2/3
    C0 = 10
    tL = 10/24
    dusk = 16/24
    lamda = 1
    H = 0.1
    b = 1
    w = 1
    K = 10
    n = 5
    S0 = 1300
    p = np.array([a,r,k,C0,tL,dusk,lamda,H,b,w,K,n,S0])
    M = model(p)
    t = M["t"].values
    C = M["C"].values
    S = M["S"].values
    F = M["F"].values
    B = M["B"].values
    plt.subplot(2,2,1)
    plt.plot(t,C)
    plt.title("C")
    plt.subplot(2,2,2)
    plt.plot(t,S)
    plt.title("S")
    plt.subplot(2,2,3)
    plt.plot(t,F)
    plt.title("F")
    plt.subplot(2,2,4)
    plt.plot(t,B)
    plt.title("B")
    plt.show()

#正規分布に従ってランダムウォークの歩幅を決める
def step(hohabas,flag):
    steps = []
    for i in range(len(flag)):
        f = flag[i]
        if f == "f":
            stp = 0
        else:
            hohaba = hohabas[i]
            stp = stats.norm.rvs(loc=0,scale=hohaba)
        steps.append(stp)
    steps = np.array(steps)
    return(steps)
#パラメータが正常な範囲内か調べ修正する
def check(p,pbs):
    m = len(p)
    for i in range(m):
        pi = p[i]
        pbi = pbs[i]
        up = pbi[1]
        low = pbi[0]
        if pi > up:
            pi  = up
        elif pi < low:
            pi = low
        p[i]=pi
    return(p)
#事前確率を計算
def prob0(p,alpha,flag):
    p0s = []
    for i in range(len(p)):
        f = flag[i]
        if f == "f":
            p0 = 1
        elif f == "n":
            par = p[i]
            mpar = alpha[2*i]
            spar = alpha[2*i+1]
            p0 = stats.norm.pdf(par,loc=mpar,scale=spar)
        else:
            par = p[i]
            uppar = alpha[2*i]
            lowpar = alpha[2*i+1]
            p0 = stats.uniform.pdf(par,loc=uppar,scale=lowpar)
        p0s.append(p0)
    p0s = np.array(p0s)
    return(p0s)
#最も近い値をリストから探す
def closest(t,tset):
    tset1 = tset - t
    tset2 = tset1**2
    min = np.min(tset2)
    t1 = tset[np.where(tset2 == min)][0]
    return(t1)
#メトロポリタンヘイスティング本体
def mcmc(df,p0,flags,pbs,alpha,stp,repeat):
    #パラメータの初期化
    p = p0
    ps = []
    #データを条件別で分割
    df_LD = df[df["Photoperiod"]=="LD"]
    df_SD = df[df["Photoperiod"]=="SD"]
    #データの時間軸？測定時点？を取得
    tset_df_LD = np.unique(df_LD["Time"].values)
    tset_df_SD = np.unique(df_SD["Time"].values)
    for i in range(repeat):
        print(i)
        ps.append(p)
        #新しいパラメータの候補
        new_p = p + step(stp,flags)
        new_p = check(new_p,pbs)
        #現在および新しいパラメータのもとでのモデルの値
        mod_LD = model(p,"LD")
        mod_SD = model(p,"SD")
        new_mod_LD = model(new_p,"LD")
        new_mod_SD = model(new_p,"SD")
        tset_mod = mod_LD["t"].values
        #採択率の初期化
        r = 1
        #LD条件でのフィットの良さ
        for t_df_LD in tset_df_LD:
            dataC = df_LD[df_LD["Time"]==t_df_LD]["Starch"].values
            dataS = df_LD[df_LD["Time"]==t_df_LD]["Sucrose"].values
            t_mod = closest(t_df_LD/24,tset_mod)
            modC = mod_LD[mod_LD["t"]==t_mod]["C"].values[0]
            modS = mod_LD[mod_LD["t"]==t_mod]["S"].values[0]
            new_modC = new_mod_LD[new_mod_LD["t"]==t_mod]["C"].values[0]
            new_modS = new_mod_LD[new_mod_LD["t"]==t_mod]["S"].values[0]
            #Starch
            for c in dataC:
                #分散
                sC = p[-2]
                #尤度
                Pc = stats.norm.pdf(c,loc=modC,scale=sC)
                #新しい分散
                new_sC = new_p[-2]
                #新しいパラメータのモデルの尤度
                new_Pc = stats.norm.pdf(c,loc=new_modC,scale=new_sC)
                #採択率の更新
                if Pc == 0:
                    rc = 1
                else:
                    rc = new_Pc/Pc
                r = r * rc
            #Sucrose
            for s in dataS:
                sS = p[-1]
                Ps = stats.norm.pdf(s,loc=modS,scale=sS)
                new_sS = new_p[-1]
                new_Ps = stats.norm.pdf(s,loc=new_modS,scale=new_sS)
                if Ps == 0:
                    rs = 1
                else:
                    rs = new_Ps/Ps
                r = r*rs
        #SD条件
        for t_df_SD in tset_df_SD:
            dataC = df_SD[df_SD["Time"]==t_df_SD]["Starch"].values
            dataS = df_SD[df_SD["Time"]==t_df_SD]["Sucrose"].values
            t_mod = closest(t_df_SD/24,tset_mod)
            modC = mod_SD[mod_SD["t"]==t_mod]["C"].values[0]
            modS = mod_SD[mod_SD["t"]==t_mod]["S"].values[0]
            new_modC = new_mod_SD[new_mod_SD["t"]==t_mod]["C"].values[0]
            new_modS = new_mod_SD[new_mod_SD["t"]==t_mod]["S"].values[0]
            #print("read data SD")
            for c in dataC:
                sC = p[-2]
                Pc = stats.norm.pdf(c,loc=modC,scale=sC)
                new_sC = new_p[-2]
                new_Pc = stats.norm.pdf(c,loc=new_modC,scale=new_sC)
                if Pc == 0:
                    rc = 1
                else:
                    rc = new_Pc/Pc
                r = r * rc
                #print("rc")
            for s in dataS:
                sS = p[-1]
                Ps = stats.norm.pdf(s,loc=modS,scale=sS)
                new_sS = new_p[-1]
                new_Ps = stats.norm.pdf(s,loc=new_modS,scale=new_sS)
                if Ps == 0:
                    rs = 1
                else:
                    rs = new_Ps/Ps
                r = r*rs
                #print("rs")
        #事前確率
        P0 = prob0(p,alpha,flags)
        new_P0 = prob0(new_p,alpha,flags)
        r0 = np.prod(new_P0/P0)
        r = r*r0
        #新パラメータを採択するか否か
        r = min(1,r)
        rand = np.random.rand()
        if rand <= r:
            p = new_p
    ps = np.array(ps).T
    return(ps)

#対数尤度
def Likelihood(df,model,params):
    simtset = model["t"].values
    datatset = np.unique(df["Time"].values)
    sC,sS = params[-2],params[-1]
    logl = 0
    for datat in datatset:
        dataCt = df[df["Time"]==datat]["Starch"].values
        dataSt = df[df["Time"]==datat]["Sucrose"].values
        simt = closest(datat/24,simtset)
        simCt = model[model["t"]==simt]["C"].values[0]
        simSt = model[model["t"]==simt]["S"].values[0]
        for c in dataCt:
            pc = stats.norm.pdf(c,loc=simCt,scale=sC)
            logpc = np.log10(pc)
            logl += logpc
        for s in dataSt:
            ps = stats.norm.pdf(s,loc=simSt,scale=sS)
            logps = np.log10(ps)
            logl += logps
    return(logl)

#結果をまとめる
def summary(df,ps,pnames,flags):
    pset = []
    pname = []
    #ｐを整理
    for i in range(len(flags)):
        flag = flags[i]
        if flag != 0:
            p = ps[i]
            pset.append(p)
            name = pnames[i]
            pname.append(name)
    #ｐのランダムウォークをプロット
    n = len(pset)
    m = int((n+1)/2)
    for i in range(n):
        plt.subplot(2,m,i+1)
        plt.plot(pset[i])
        plt.title(pname[i])
    plt.show()
    #バーンインを設定
    barnin = int(input("Please set burnin >> "))
    #ヒストグラム
    for i in range(n):
        plt.subplot(2,m,i+1)
        p = pset[i]
        l = len(p[barnin:])
        N = int(l**0.5)
        plt.hist(p[barnin:],bins = N)
        plt.title(pname[i])
    plt.show()
    #統計量の算出
    pmeans = np.mean(pset[:],axis=1)
    pstds = np.std(pset[:],axis=1)
    percentiles975 = [np.percentile(pset[i],97.5) for i in range(n)]
    percentiles225 = [np.percentile(pset[i],2.25) for i in range(n)]
    #まとめ
    sum = pd.DataFrame(index=pname)
    sum["mean"] = pmeans
    sum["std"] = pstds
    sum["97.5%"] = percentiles975
    sum["2.25%"] = percentiles225
    print(sum)
    #データとの比較
    ps_all = ps.T[barnin:]
    print(len(ps_all))
    pmeans_all = np.mean(ps_all,axis=0)
    print(len(pmeans_all))
    sim_LD = model(pmeans_all,"LD")
    sim_SD = model(pmeans_all,"SD")
    #尤度
    like_LD = Likelihood(df,sim_LD,pmeans_all)
    like_SD = Likelihood(df,sim_SD,pmeans_all)
    like = like_LD + like_SD
    print("Likelihood LD:",like_LD," SD:",like_SD," total:",like)
    simt_LD = sim_LD["t"].values*24
    simC_LD = sim_LD["C"].values
    simS_LD = sim_LD["S"].values
    simF_LD = sim_LD["F"].values*24
    simB_LD = sim_LD["B"].values
    data_LD = df[df["Photoperiod"]=="LD"]
    datat_LD = data_LD["Time"].values
    dataC_LD = data_LD["Starch"].values
    dataS_LD = data_LD["Sucrose"].values
    dataB_LD = data_LD["Beta"].values
    plt.subplot(2,4,1)
    plt.plot(simt_LD,simC_LD)
    plt.scatter(datat_LD,dataC_LD)
    plt.title("Starch LD")
    plt.subplot(2,4,2)
    plt.plot(simt_LD,simS_LD)
    plt.scatter(datat_LD,dataS_LD)
    plt.title("Sucrose LD")
    plt.subplot(2,4,3)
    plt.plot(simt_LD,simF_LD)
    plt.title("Phi LD")
    plt.subplot(2,4,4)
    plt.plot(simt_LD,simB_LD)
    plt.scatter(datat_LD,dataB_LD)
    plt.title("Beta LD")

    simt_SD = sim_SD["t"].values*24
    simC_SD = sim_SD["C"].values
    simS_SD = sim_SD["S"].values
    simF_SD = sim_SD["F"].values*24
    simB_SD = sim_SD["B"].values
    data_SD = df[df["Photoperiod"]=="SD"]
    datat_SD = data_SD["Time"].values
    dataC_SD = data_SD["Starch"].values
    dataS_SD = data_SD["Sucrose"].values
    dataB_SD = data_SD["Beta"].values
    plt.subplot(2,4,5)
    plt.plot(simt_SD,simC_SD)
    plt.scatter(datat_SD,dataC_SD)
    plt.title("Starch SD")
    plt.subplot(2,4,6)
    plt.plot(simt_SD,simS_SD)
    plt.scatter(datat_SD,dataS_SD)
    plt.title("Sucrose SD")
    plt.subplot(2,4,7)
    plt.plot(simt_SD,simF_SD)
    plt.title("Phi SD")
    plt.subplot(2,4,8)
    plt.plot(simt_SD,simB_SD)
    plt.scatter(datat_SD,dataB_SD)
    plt.title("Beta SD")
    plt.show()
    return(sum)

#実行関数
#パラメータの初期設定等はこの関数の内部でできます。
def sim_wt():
    ## パラメータ
    pnames = ["a","r","k","C0_LD","C0_LD","tL","dusk_LD","dusk_SD","lamda","H","b","w","K","n","S0_LD","S0_SD","sC","sS"]
    # fixed
    a = 144
    k = 2/3
    C0_LD = 9.996793
    C0_SD = 5.779395
    tL = 10/24
    dusk_LD = 16/24
    dusk_SD = 8/24
    lamda = 1 #(wt)
    b = 1
    w = 1
    S0_LD = 1037
    S0_SD = 549.300137
    # free
    r = 0.8
    H = 0.5
    K = 10
    n = 20
    sC = 5
    sS = 500
    p0 = np.array([a,r,k,C0_LD,C0_SD,tL,dusk_LD,dusk_SD,lamda,H,b,w,K,n,S0_LD,S0_SD,sC,sS])
    # flags:パラメータを変化させるかどうかを指定
    ##　ｆ　＝　固定、ｎ＝事前分布が正規分布、uni＝事前確率が一様分布
    flags = ["f","n","f","f","f","f","f","f","f","n","f","f","n","n","f","f","n","n"]
    #         a   r   k   cl cs  tL   dl  ds  l   H   b   w   K   n   sl  ss  sc  ss
    # alpha:事前確率の平均と分散（正規分布）または下限と上限（一様分布）
    ma,sa = 144,50
    mr,sr = 0.5,0.4
    mk,sk = 2/3,0.1
    mC0,sC0 = 8.127,3.20
    mtL,stL = 10/24,0.4
    mdL,sdL = 18/24,0.1
    mdS,sdS = 8/24,0.1
    ml,sl = 1,0.5
    mH,sH = 0.5,0.4
    mb,sb = 1,0.5
    mw,sw = 1,0.5
    mK,sK = 20,5
    mn,sn = 6,1
    mS0,sS0 = 820,388.6
    msC,ssC = 15,3.3
    msS,ssS = 494.65,43.86
    alpha = np.array([ma,sa,mr,sr,mk,sk,mC0,sC0,mC0,sC0,mtL,stL,mdL,sdL,mdS,sdS,ml,sl,mH,sH,mb,sb,mw,sw,mK,sK,mn,sn,mS0,sS0,mS0,sS0,msC,ssC,msS,ssS])
    #pbsパラメータの範囲（pb+parameter）
    M = 100000
    pba = [0,M]
    pbk = [0,1]
    pbC0 = [0,M]
    pbtL = [0,1]
    pbdusk_LD = [0,1]
    pbdusk_SD = [0,1]
    pblamda = [0,M]
    pbb = [0,M]
    pbw = [0,M]
    pbS0 = [0,M]
    pbr = [0,1]
    pbH = [0,1]
    pbK = [0,M]
    pbn = [0,M]
    pbsC = [0,M]
    pbsS = [0,M]
    pbs = np.array([pba,pbr,pbk,pbC0,pbC0,pbtL,pbdusk_LD,pbdusk_SD,pblamda,pbH,pbb,pbw,pbK,pbn,pbS0,pbS0,pbsC,pbsS])
    #stepのスケール
    stpa = 0.01
    stpr = 0.01
    stpk = 0.01
    stpC0 = 0.02
    stptL = 0.01
    stpdusk_LD = 0
    stpdusk_SD = 0
    stplamda = 0.01
    stpH = 0.02
    stpb = 0.01
    stpw = 0.01
    stpK = 0.1
    stpn = 0.1
    stpS0 = 0.1
    stpsC = 0.05
    stpsS = 1
    stp = np.array([stpa,stpr,stpk,stpC0,stpC0,stptL,stpdusk_LD,stpdusk_SD,stplamda,stpH,stpb,stpw,stpK,stpn,stpS0,stpS0,stpsC,stpsS])
    #mcmc
    repeat = 20000
    df = pd.read_csv("wt_fitting_data.csv")
    ps = mcmc(df,p0,flags,pbs,alpha,stp,repeat)
    result = summary(df,ps,pnames,flags)


sim_wt()
