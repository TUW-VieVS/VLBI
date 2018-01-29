% ************************************************************************
%   Description:
%   function to exclude  models of right ascension and declination as sourcewise 
%
%   Reference: 
%
%   Input:										
%       'ns'        (1,1)                number of sources
%       'tso'        structure array     estimation intervals of source coor. for each source
%       'opt'        structure array     (for info. /DOC/opt.doc)
%       'nso'        structure array     the number of source coordinate offsets for each source
%       'sumso'      structure array     % total vector of source coor. estimates after eleminating
                                         % non-observed sources
%       'A'          structure array     design matrix of the real-observation equations
%       'H'          structure array     design matrix of the pseudo-observation equations (constraints)
%       'Ph'         structure array     weight matrix of the pseudo-observation equations 
%       'och'        structure array     o-c vector of constraints (zero vector)
%
%   Output:
%       'nso'        structure array     the number of source coordinate offsets for each source
%       'tso'        structure array     estimation intervals of source coor. for each source
%       'A'          structure array     design matrix of the real-observation equations
%       'H'          structure array     design matrix of the pseudo-observation equations (constraints)
%       'Ph'         structure array     weight matrix of the pseudo-observation equations 
%       'och'        structure array     o-c vector of constraints (zero vector)
%
%   External calls: 	
%   
%   Coded for VieVS: 
%   26 Aug 2009 by Kamil Teke
%
%   Revision: 
%   06 Dec 2009 by Kamil Teke: header added
%   18 Jun 2014 by Hana Krasna: sources can be estimated in vie_lsm with
%   NNR condition
% ************************************************************************
function [nso,tso,A,H,Ph,och] = delsource(ns,tso,opt,nso,sumso,A,H,Ph,och)

count = 0;
del = 0;
vecm = 0:1:ns;
if opt.est_sourceNNR==0 
    for isou = 1 : ns 
        if opt.source(ns+1-isou).rade_inc == 0

            count = count + 1;
            del(count) = ns+1-isou;

            summ_rade = sumso.sources-vecm;
            A(11).sm(:,sumso.sources(del(count))+1:sumso.sources(del(count)+1)) = [];
            A(12).sm(:,sumso.sources(del(count))+1:sumso.sources(del(count)+1)) = [];

            H(11).sm(:,sumso.sources(del(count))+1:sumso.sources(del(count)+1)) = [];
            H(11).sm(summ_rade(del(count))+1:summ_rade(del(count)+1),:) = [];

            H(12).sm = [];
            H(12).sm = H(11).sm; 

            Ph(11).sm(:,summ_rade(del(count))+1:summ_rade(del(count)+1)) = []; 
            Ph(11).sm(summ_rade(del(count))+1:summ_rade(del(count)+1),:) = [];

            Ph(12).sm = []; 
            Ph(12).sm = Ph(11).sm; 

            if opt.constr_sou == 1
                och(12).sv = []; 
                och(11).sv(summ_rade(del(count))+1:summ_rade(del(count)+1)) = [];
                och(12).sv = och(11).sv; 
            end

            nso(:,del(count)).sources = 0;
            tso(del(count)) = [];
        end
    end
end