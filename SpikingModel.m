%%%%%
%
% Implementation of network of spiking neurons
%
%%%%%

clear;

% Duration of simulation
TmTEnd = 2; % Seconds

% Time step
EvltnTmStp = 0.0001;


function [Thrshlds, RstPtntls, RfrctryPrd, EPSPDrtn, IPSPDrtn, EPSPAmpltd, IPSPAmpltd, InptPrbblty, IPSPInptPrbblty, NmbrNrnsVctr, NmbrExcttryPostSynptc, ExcttryPstSynptcClls, NmbrInhbtryPostSynptc, InhbtryPstSynptcClls, FrmPyrmdlCllsToExcttryIntrnrns, FrmExcttryIntrnrnsTPyrmdlClls, FrmPyrmdlCllsToInhbtryIntrnrns, FrmInhbtryIntrnrnsTPyrmdlClls] = GnrtNtwrkPrmtrs(EvltnTmStp)

  % Network parameters, as defined in paper text

  EPSPDrtn = 100; % Reciprocal of EPSP time constant 
  IPSPDrtn = 50; % Reciprocal of IPSP time constant 

  EPSPAmpltd = 0.00325; % EPSP amplitude
  IPSPAmpltd = -0.022; % IPSP amplitude

  ThrshldFdFrwrdPyrmdlClls = 0.01352; % Action potential threshold for feedforward pyramidal cells
  ThrshldExcttryIntrnrns = 1000000; % Action potential threshold for excitatory interneurons
  ThrshldInhbtryIntrnrns = 0.004956; % Action potential threshold for inhibitory interneuons

  % Note: Number of neurons, connection probability, and number of synpases are linked.

  NmbrPyrmdlNrns = 750; % Number of feedforward pyramidal cells
  NmbrLclExcttryNrns = 2; % Number of excitatory interneurons
  NmbrLclInhbtryNrns = 750; % Number of inhibitory interneuons

  %% Number of connections recevied per postsynaptic cell:
  FrmPyrmdlCllsToExcttryIntrnrns = 0; % From feedforward pyramidal cells to excitatory interneurons
  FrmExcttryIntrnrnsTPyrmdlClls = 0; % From excitatory interneurons to feedforward pyramidal cells
  FrmPyrmdlCllsToInhbtryIntrnrns = 75; % From feedforward pyramidal cells to inhibitory interneurons
  FrmInhbtryIntrnrnsTPyrmdlClls = 75; % From nhibitory interneurons to feedforward pyramidal cells

  NmbrNrnsVctr(1) = NmbrPyrmdlNrns;
  NmbrNrnsVctr(2) = NmbrLclExcttryNrns;
  NmbrNrnsVctr(3) = NmbrLclInhbtryNrns;

  TtlNmbrNrns = NmbrPyrmdlNrns + NmbrLclExcttryNrns + NmbrLclInhbtryNrns;

  NrnCV = 1;

  for n=1:NmbrPyrmdlNrns

      Thrshlds(NrnCV, 1) = ThrshldFdFrwrdPyrmdlClls;

      NrnCV = NrnCV + 1;

  end

  for n=1:NmbrLclExcttryNrns

      Thrshlds(NrnCV, 1) = ThrshldExcttryIntrnrns;

      NrnCV = NrnCV + 1;

  end

  for n=1:NmbrLclInhbtryNrns

      Thrshlds(NrnCV, 1) = ThrshldInhbtryIntrnrns;

      NrnCV = NrnCV + 1;

  end

  RstPtntls = zeros(TtlNmbrNrns, 1);
  RstPtntls = RstPtntls - 0;

  RfrctryPrd = zeros(TtlNmbrNrns, 1);
  RfrctryPrd = RfrctryPrd + 0.001;

  NrnCV = 1;

  for n=1:NmbrPyrmdlNrns

      InptPrbblty(NrnCV) = 220*EvltnTmStp;
      IPSPInptPrbblty(NrnCV) = 0;

      NrnCV = NrnCV + 1;

  end

  for n=1:NmbrLclExcttryNrns

      InptPrbblty(NrnCV) = 0;
      IPSPInptPrbblty(NrnCV) = 0;

      NrnCV = NrnCV + 1;

  end

  for n=1:NmbrLclInhbtryNrns

      InptPrbblty(NrnCV) = 0;
      IPSPInptPrbblty(NrnCV) = 0;

      NrnCV = NrnCV + 1;

  end

  CnnctnWghts = zeros(TtlNmbrNrns);

  PyrmdlNrnsStrt = 1;
  LclExcttryNrnsStrt = PyrmdlNrnsStrt + NmbrPyrmdlNrns;
  LclInhbtryNrnsStrt = LclExcttryNrnsStrt + NmbrLclExcttryNrns;

  for n=1:NmbrLclExcttryNrns
    CntNrns = FrmPyrmdlCllsToExcttryIntrnrns;
    while (CntNrns>0)
      NrnNmbr = ceil(rand*NmbrPyrmdlNrns);
      if (CnnctnWghts(LclExcttryNrnsStrt + n - 1, PyrmdlNrnsStrt + NrnNmbr - 1) == 0)
        CnnctnWghts(LclExcttryNrnsStrt + n - 1, PyrmdlNrnsStrt + NrnNmbr - 1) = 1;
        CntNrns = CntNrns - 1;
       end
    end
  end

  for n=1:NmbrPyrmdlNrns
    CntNrns = FrmExcttryIntrnrnsTPyrmdlClls;
    while (CntNrns>0)
      NrnNmbr = ceil(rand*NmbrLclExcttryNrns);
      if (CnnctnWghts(PyrmdlNrnsStrt + n - 1, LclExcttryNrnsStrt + NrnNmbr - 1) == 0)
        CnnctnWghts(PyrmdlNrnsStrt + n - 1, LclExcttryNrnsStrt + NrnNmbr - 1) = 1;
        CntNrns = CntNrns - 1;
       end
    end
  end

  for n=1:NmbrLclInhbtryNrns
    CntNrns = FrmPyrmdlCllsToInhbtryIntrnrns;
    while (CntNrns>0)
      NrnNmbr = ceil(rand*NmbrPyrmdlNrns);
      if (CnnctnWghts(LclInhbtryNrnsStrt + n - 1, PyrmdlNrnsStrt + NrnNmbr - 1) == 0)
        CnnctnWghts(LclInhbtryNrnsStrt + n - 1, PyrmdlNrnsStrt + NrnNmbr - 1) = 1;
        CntNrns = CntNrns - 1;
       end
    end
  end

  for n=1:NmbrPyrmdlNrns
    CntNrns = FrmInhbtryIntrnrnsTPyrmdlClls;
    while (CntNrns>0)
      NrnNmbr = ceil(rand*NmbrLclInhbtryNrns);
      if (CnnctnWghts(PyrmdlNrnsStrt + n - 1, LclInhbtryNrnsStrt + NrnNmbr - 1) == 0)
        CnnctnWghts(PyrmdlNrnsStrt + n - 1, LclInhbtryNrnsStrt + NrnNmbr - 1) = -1;
        CntNrns = CntNrns - 1;
       end
    end
  end

  NmbrExcttryPostSynptc = zeros(sum(NmbrNrnsVctr), 1);
  ExcttryPstSynptcClls = zeros(sum(NmbrNrnsVctr), 1);
  for n=1:sum(NmbrNrnsVctr)
    for m=1:sum(NmbrNrnsVctr)
      if (CnnctnWghts(m,n) == 1)
        NmbrExcttryPostSynptc(n) = NmbrExcttryPostSynptc(n) + 1;
        ExcttryPstSynptcClls(n, NmbrExcttryPostSynptc(n)) = m;
      end
    end
  end

  NmbrInhbtryPostSynptc = zeros(sum(NmbrNrnsVctr), 1);
  InhbtryPstSynptcClls = zeros(sum(NmbrNrnsVctr), 1);
  for n=1:sum(NmbrNrnsVctr)
    for m=1:sum(NmbrNrnsVctr)
      if (CnnctnWghts(m,n) == -1)
        NmbrInhbtryPostSynptc(n) = NmbrInhbtryPostSynptc(n) + 1;
        InhbtryPstSynptcClls(n, NmbrInhbtryPostSynptc(n)) = m;
      end
    end
  end

end


function [Tm, LFP, SpkCnt, ExmplTrcs] = TmEvltn(Thrshlds, RstPtntls, RfrctryPrd, EPSPDrtn, IPSPDrtn, EPSPAmpltd, IPSPAmpltd, InptPrbblty, IPSPInptPrbblty, NmbrNrnsVctr, EvltnTmStp, TmTEnd, NmbrExcttryPostSynptc, ExcttryPstSynptcClls, NmbrInhbtryPostSynptc, InhbtryPstSynptcClls)

  [NmbrNrns, Vrbls] = size(Thrshlds);

  TmLstSpk = zeros(NmbrNrns, 1);
  TmLstSpk = TmLstSpk - 100;

  SpkCnt = zeros(NmbrNrns, 1);

  InptStps = 2;

  t = 0;
  V = RstPtntls;

  EPSPs(1:NmbrNrns,1) = 0;
  IPSPs(1:NmbrNrnsVctr(1),1) = 0;

  dEPSPs(1:NmbrNrns,1) = 0;
  dIPSPs(1:NmbrNrnsVctr(1),1) = 0;

  Inpt = zeros(NmbrNrns,InptStps);

  InhbtryInpt = zeros(NmbrNrnsVctr(1), InptStps);
  ExcttryInpt = zeros(NmbrNrns, InptStps);

  for n=1:NmbrNrnsVctr(1)
    for m=1:InptStps
      if (rand < InptPrbblty(n))
        ExcttryInpt(n,m) = 1/EvltnTmStp;
      else
        ExcttryInpt(n,m) = 0;
      end
    end
  end

  SvIndx = 1;

  while (t <= TmTEnd)

      for n=1:NmbrNrnsVctr(1)
          if (TmLstSpk(n) + RfrctryPrd(n) >= t)
            dEPSPs(n,1) = 0;
            dIPSPs(n,1) = 0;
            EPSPs(n,1) = 0;
            IPSPs(n,1) = 0;
            V(n) = V(n) - EvltnTmStp*(Thrshlds(n) - RstPtntls(n))/RfrctryPrd(n);
          else
            dV = 0;
            dEPSPs(n,1) = dEPSPs(n,1) + EvltnTmStp*(EPSPAmpltd*EPSPDrtn*ExcttryInpt(n,1) - 2*EPSPDrtn*dEPSPs(n,1) - EPSPDrtn^2*EPSPs(n,1));
            EPSPs(n,1) = EPSPs(n,1) + EvltnTmStp * dEPSPs(n,1);
            dIPSPs(n,1) = dIPSPs(n,1) + EvltnTmStp*(IPSPAmpltd*IPSPDrtn*InhbtryInpt(n,1) - 2*IPSPDrtn*dIPSPs(n,1) - IPSPDrtn^2*IPSPs(n,1));
            IPSPs(n,1) = IPSPs(n,1) + EvltnTmStp * dIPSPs(n,1);
            V(n) = RstPtntls(n) + EPSPs(n,1) + IPSPs(n,1);
          end
      end

      for n=NmbrNrnsVctr(1)+1:NmbrNrns
          if (TmLstSpk(n) + RfrctryPrd(n) >= t)
            dEPSPs(n,1) = 0;
            EPSPs(n,1) = 0;
            V(n) = V(n) - EvltnTmStp*(Thrshlds(n) - RstPtntls(n))/RfrctryPrd(n);
          else
            dV = 0;
            dEPSPs(n,1) = dEPSPs(n,1) + EvltnTmStp*(EPSPAmpltd*EPSPDrtn*ExcttryInpt(n,1) - 2*EPSPDrtn*dEPSPs(n,1) - EPSPDrtn^2*EPSPs(n,1));
            EPSPs(n,1) = EPSPs(n,1) + EvltnTmStp * dEPSPs(n,1);
            V(n) = RstPtntls(n) + EPSPs(n,1);
          end
      end


      for n=1:NmbrNrnsVctr(1)+NmbrNrnsVctr(2)
          if (V(n)>=Thrshlds(n) && t > TmLstSpk(n) + RfrctryPrd(n))
              for m=1:NmbrExcttryPostSynptc(n)
                  ExcttryInpt(ExcttryPstSynptcClls(n,m), 2) = ExcttryInpt(ExcttryPstSynptcClls(n,m), 2) + 1/EvltnTmStp;
              end

              SpkCnt(n) = SpkCnt(n) + 1;

              TmLstSpk(n) = t;
          end
      end

      for n=NmbrNrnsVctr(1)+NmbrNrnsVctr(2)+1:NmbrNrns
          if (V(n)>=Thrshlds(n) && t > TmLstSpk(n) + RfrctryPrd(n))
              for m=1:NmbrInhbtryPostSynptc(n)
                  InhbtryInpt(InhbtryPstSynptcClls(n,m), 2) = InhbtryInpt(InhbtryPstSynptcClls(n,m), 2) + 1/EvltnTmStp;
              end

              SpkCnt(n) = SpkCnt(n) + 1;

              TmLstSpk(n) = t;
          end
      end


      Tm(SvIndx) = t;
      LFP(SvIndx) = mean(V(1:NmbrNrnsVctr(1)));

      ExmplTrcs(1,SvIndx) = V(1);
      ExmplTrcs(2,SvIndx) = V(2);

      ExmplTrcs(3,SvIndx) = V(NmbrNrnsVctr(1)+1);
      ExmplTrcs(4,SvIndx) = V(NmbrNrnsVctr(1)+2);

      ExmplTrcs(5,SvIndx) = V(NmbrNrnsVctr(1)+NmbrNrnsVctr(2)+1);
      ExmplTrcs(6,SvIndx) = V(NmbrNrnsVctr(1)+NmbrNrnsVctr(2)+2);

      SvIndx = SvIndx + 1;


      ExcttryInpt = circshift(ExcttryInpt, -1, 2);
      InhbtryInpt = circshift(InhbtryInpt, -1, 2);

      ExcttryInpt(:,InptStps) = 0;
      InhbtryInpt(:,InptStps) = 0;

      for n=1:NmbrNrnsVctr(1)
        if (rand < InptPrbblty(n))
          ExcttryInpt(n,InptStps) = 1/EvltnTmStp;
        end
      end

      t = t + EvltnTmStp;

  end

end

function [OutputTF] = NrmlsdTrnsfrFnctn(Prmtrs, fInptFT)

  a = Prmtrs(1,1);
  A = Prmtrs(2,1);

  b = Prmtrs(3,1);
  B = Prmtrs(4,1);

  r = Prmtrs(5,1);
  e0 = Prmtrs(6,1);
  v0 = Prmtrs(7,1);

  C1 = Prmtrs(8,1);
  C2 = Prmtrs(9,1);
  C3 = Prmtrs(10,1);
  C4 = Prmtrs(11,1);

  gmm = (2*e0*r*exp(r*v0))/((1+exp(r*v0))^2);

  [vrbl0 NmbrPnts] = size(fInptFT);

  Frqncs = fInptFT;

  %%%
  %
  % Note the factor of 2pi introduced here to equate with FFT
  %
  %%%

  for n=1:NmbrPnts

      s = 1i * 2 * pi * Frqncs(n);

      Trm1 = ((s+a)^2)/(A*a);
      Trm2 = (a*A*C1*C2*gmm^2)/((s+a)^2);
      Trm3 = (b*B*C3*C4*gmm^2)/((s+b)^2);

      TF(n) = 1/(Trm1 - Trm2 + Trm3);
      OutputTF(n) = TF(n);
      OutputTF(n) = (abs(OutputTF(n)));
      fTF(n) = fInptFT(n);

  end

  OutputTF = OutputTF/mean(OutputTF);

end

[Thrshlds, RstPtntls, RfrctryPrd, EPSPDrtn, IPSPDrtn, EPSPAmpltd, IPSPAmpltd, InptPrbblty, IPSPInptPrbblty, NmbrNrns, NmbrExcttryPostSynptc, ExcttryPstSynptcClls, NmbrInhbtryPostSynptc, InhbtryPstSynptcClls, C1, C2, C3, C4] = GnrtNtwrkPrmtrs(EvltnTmStp);

[Tm, LFP, SpkCnt, ExmplTrcs] = TmEvltn(Thrshlds, RstPtntls, RfrctryPrd, EPSPDrtn, IPSPDrtn, EPSPAmpltd, IPSPAmpltd, InptPrbblty, IPSPInptPrbblty, NmbrNrns, EvltnTmStp, TmTEnd, NmbrExcttryPostSynptc, ExcttryPstSynptcClls, NmbrInhbtryPostSynptc, InhbtryPstSynptcClls);

[vrbl0, NmbrPnts] = size(LFP);

FT(1,:) = fft(LFP);
FT(1,:) = fftshift(FT(1,:));
fFll = (-(NmbrPnts-1)/2:(NmbrPnts-1)/2)*1/EvltnTmStp*1/NmbrPnts;
SvFT(1,:) = abs(FT(1,:)/NmbrPnts);
[vrbl0, IndxMn] = min(abs(fFll-1));
[vrbl0, IndxMx] = min(abs(fFll-100));
Svf = fFll(IndxMn:IndxMx);
SvFT1 = SvFT(1,IndxMn:IndxMx);

Prmtrs(1,1) = EPSPDrtn;
Prmtrs(2,1) = EPSPAmpltd;
Prmtrs(3,1) = IPSPDrtn;
Prmtrs(4,1) = -1*IPSPAmpltd;
Prmtrs(5,1) = 560;
Prmtrs(6,1) = 2.5;
Prmtrs(7,1) = 0.006;
Prmtrs(8,1) = C1;
Prmtrs(9,1) = C2;
Prmtrs(10,1) = C3;
Prmtrs(11,1) = C4;

TF = NrmlsdTrnsfrFnctn(Prmtrs, Svf);

FntSz = 25;

figure(1);
subplot(3,2,1);
hold on;
plot(Tm, 1000*ExmplTrcs(1,:), 'k', 'linewidth', 2);
plot(Tm, 1000*ExmplTrcs(2,:), 'k', 'linewidth', 2);
set(gca, 'fontsize', FntSz);
title('Example traces: Feedforward neurons', 'fontsize', FntSz);
xlabel('Time (s)', 'fontsize', FntSz);
ylabel('Potential (mV)', 'fontsize', FntSz);

subplot(3,2,3);
hold on;
plot(Tm, 1000*ExmplTrcs(3,:), 'k', 'linewidth', 2);
plot(Tm, 1000*ExmplTrcs(4,:), 'k', 'linewidth', 2);
set(gca, 'fontsize', FntSz);
title('Example traces: Excitatory interneurons', 'fontsize', FntSz);
xlabel('Time (s)', 'fontsize', FntSz);
ylabel('Potential (mV)', 'fontsize', FntSz);

subplot(3,2,5);
hold on;
plot(Tm, 1000*ExmplTrcs(5,:), 'k', 'linewidth', 2);
plot(Tm, 1000*ExmplTrcs(6,:), 'k', 'linewidth', 2);
set(gca, 'fontsize', FntSz);
title('Example traces: Inhibitory interneurons', 'fontsize', FntSz);
xlabel('Time (s)', 'fontsize', FntSz);
ylabel('Potential (mV)', 'fontsize', FntSz);

subplot(3,2,2);
plot(Tm, 1000*LFP, 'k', 'linewidth', 2);
set(gca, 'fontsize', FntSz);
title('Effective local field potential', 'fontsize', FntSz);
xlabel('Time (s)', 'fontsize', FntSz);
ylabel('Potential (mV)', 'fontsize', FntSz);

subplot(3,2,4);
plot(Svf, SvFT1/mean(SvFT1), 'k', 'linewidth', 2);
set(gca, 'fontsize', FntSz);
title('Fourier transform of local field potential', 'fontsize', FntSz);
xlabel('Frequency (Hz)', 'fontsize', FntSz);
ylabel('Amplitude', 'fontsize', FntSz);

subplot(3,2,6);
plot(Svf, TF/mean(TF), 'k', 'linewidth', 2);
set(gca, 'fontsize', FntSz);
title('Equivalent transfer function', 'fontsize', FntSz);
xlabel('Frequency (Hz)', 'fontsize', FntSz);
ylabel('Amplitude', 'fontsize', FntSz);
