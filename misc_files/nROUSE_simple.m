function nROUSE_simple

% assignments for visual layer

    VPR=1;   % visual prime
    VTR=2;   % visual target
    VMK=3;   % visual mask
    VTRC=4;  % visual target choice
    VFLC=5;  % visual foil choice

% assignments for orthographic and semantic layers

    TARG=1;  % target
    FOIL=2;  % foil
    
% weight matrices

    VisOrth=[0,0;  % from VPR
             1,0;  % from VTR
             0,0;  % from VMK
             1,0;  % from VTRC
             0,1]; % from VFLC
         
    OrthSem=eye(2); % identity matrix for weights
    SemOrth=eye(2); % same identity matrix for feedback
    
% experimental procedures

    TarDur=50;
    MaskDur=500-TarDur;
    ChoiceDur=500;
    durations=[17,50,150,400,2000];
    
% parameters

    F=.25;          % semantic to orthographic feedback scalar
    N=0.0302;    % noise constant
    L=.15;          % constant leak current
    D=.324;         % synaptic depletion rate
    R=.022;         % recovery rate
    I=0.9844;       % inhibition constant
    T=.15;          % activation threshold
    S=[0.0294,0.0609,.015];    % integration time constants at each level
    
    
    for cd=1:2      % cd is a step through index for conditions    
        if cd==1         % target primed
            VisOrth(VPR,:)=[2,0];              % set to 2 because there are two visual copies on the screen
        elseif cd==2      % foil primed
            VisOrth(VPR,:)=[0,2];
        end

        [accs(:,cd), Latency]=simulate;  % run all prime durations and return accuracy and latency
        targ_lat(:,cd)=Latency(:,TARG);
        foil_lat(:,cd)=Latency(:,FOIL);
    end
    subplot(3,1,1);
    accs
    plot(log(durations),accs);
    subplot(3,1,2);
    targ_lat
    plot(log(durations),targ_lat);
    subplot(3,1,3);
    foil_lat
    plot(log(durations),foil_lat);

    
    function [acc, Latency]=simulate
        
        Latency=zeros(size(durations,2),2);     % identification latencies for each choice
        for pd=1:size(durations,2)      % pd is a step through index for prime durations
            
            PrimeDur=durations(pd);
            SOA=PrimeDur+TarDur+MaskDur;        % time when choices are presented
            
            mem_vis=zeros(1,5);     % initialize neural variables
            amp_vis=ones(1,5);
            out_vis=zeros(1,5);
            mem_orth=zeros(1,2);
            amp_orth=ones(1,2);
            out_orth=zeros(1,2);
            mem_sem=zeros(1,2);
            amp_sem=ones(1,2);
            out_sem=zeros(1,2);
            old_sem=zeros(1,2);     % needed to check for peak output


            for t=1:PrimeDur+TarDur+MaskDur+ChoiceDur
                
                % udpate visual layer
                if t==1                         % present prime
                    inp_vis=zeros(1,5);
                    inp_vis(VPR)=1;
                elseif t==PrimeDur+1           % present target
                    inp_vis=zeros(1,5);
                    inp_vis(VTR)=1;
                elseif t==PrimeDur+TarDur+1    % present mask
                    inp_vis=zeros(1,5);
                    inp_vis(VMK)=1;
                elseif t==SOA+1                 % present choices
                    inp_vis=zeros(1,5);
                    inp_vis(VTRC)=1;
                    inp_vis(VFLC)=1;
                end
                [new_mem_vis,new_amp_vis,out_vis]=update(mem_vis,amp_vis,inp_vis,1);
                
                % update orthographic layer
                inp_orth=out_vis*VisOrth;
                inp_orth=inp_orth+F.*out_sem*SemOrth;
                [new_mem_orth,new_amp_orth,out_orth]=update(mem_orth,amp_orth,inp_orth,2);

                % update semantic layer
                inp_sem=out_orth*OrthSem;
                [new_mem_sem,new_amp_sem,out_sem]=update(mem_sem,amp_sem,inp_sem,3);       

                % perceptual decision process
                if t>SOA+50    % the +50 gives things a chance to get going before peak activation is checked
                    for tf=1:2  % step through index for target and foil
                        if out_sem(tf)<old_sem(tf) && Latency(pd,tf)==0    % check for peak activation
                            Latency(pd,tf)=t-SOA;
                        end
                    end
                    old_sem=out_sem;
                end
                
                % swap new variables for old variables
                
                mem_vis=new_mem_vis;
                amp_vis=new_amp_vis;
                mem_orth=new_mem_orth;
                amp_orth=new_amp_orth;
                mem_sem=new_mem_sem;
                amp_sem=new_amp_sem;

            end
            
            % calculate accuracy
            mean_diff=Latency(pd,FOIL)-Latency(pd,TARG); % average difference between target and foil latency
            var_diff=sum(exp(N.*Latency(pd,:)));         % variance of difference between target and foil latency
            
            acc(pd,1)=1-normcdf(0,mean_diff,var_diff.^.5);
            
            if Latency(pd,TARG)==0 && Latency(pd,FOIL)>0      % target never launched
                acc(pd,1)=0;
            elseif Latency(pd,TARG)>0 && Latency(pd,FOIL)==0  % foil never launched
                acc(pd,1)=1;
            elseif Latency(pd,TARG)==0 && Latency(pd,FOIL)==0 % neither launched
                acc(pd,1)=.5;
            end
        end
    end


    function [new_mem,new_amp,old_out]=update(old_mem,old_amp,inp,level)

        old_out=subplus(old_mem-T).*old_amp;    % output is above threshold activation times available synaptic resources

        if level==1     % visual inhibition                 
            inhibit([1:3])=sum(old_out([1:3]));  % prime, target, and mask mutually inhibit each other (centrally presented)
            inhibit([4:5])=old_out([4:5]);      % choice words only self inhibit (unique screen location)
        else
            inhibit=sum(old_out);               % for orthography and semantics, everything inhibits everything
        end
                    % update membrane potential and synaptic resources
                    
        new_mem=old_mem + (S(level) .*  (  ( inp .* (1 - old_mem) )  -  (L.*old_mem)  -  (I.*inhibit.*old_mem)  ));
        new_amp=old_amp + (S(level) .*  (  (  R  .* (1 - old_amp) )  -  (D.*old_out))  );

        new_mem(new_mem<0)=0;       % keep things within bounds
        new_mem(new_mem>1)=1;
        new_amp(new_amp<0)=0;
        new_amp(new_amp>1)=1;
    end
end
