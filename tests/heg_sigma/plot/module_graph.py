#!/usr/bin/env python
__author__ = "Tommaso Chiarotti"
__license__ = "see LICENSE.txt file."
__version__ = "0"
__email__ = "tommaso.chiarotti@epfl.ch"
"""
This module contains function to plot in a systematic way the outputs of heg_sigma.x
"""

import numpy as np
import matplotlib.pyplot as plt

class plot:
    """
    this class is the plottng class 
    """
    def __init__(self,dirfile,listiter,listdelta,listpps,vstart,vend,vstep,xstart,xend,xstep):
        """
        this function initialize the plot class
        """
        self.dirfile   = dirfile
        self.listiter  = listiter 
        self.listdelta = listdelta
        self.listpps   = listpps
        self.vstart    = vstart
        self.vend      = vend
        self.vstep     = vstep
        self.xstart    = xstart
        self.xend      = xend
        self.xstep     = xstep

    def sigma_and_greenf(self,xylabel,xylim,fileout='conv_tests',classfile=['sigma_plot','greenfunc_plot']):
        """
        this functions plots the self-energy and the greensfunction
        """
        for k in self.listpps:
            #
            for iiter in self.listiter:
                #
                for nq in range(self.vstart,self.vend,self.vstep):
                    #
                    fig = plt.figure()
                    #
                    ax = []
                    for iax in range(0,4):
                        ax.append(fig.add_subplot(4,1,iax+1))
                    #
                    ax[0].set_ylabel('Re(Sigma)')
                    #ax[0].set_xlim(-10,10)
                    ax[1].set_ylabel('Im(Sigma)')
                    #ax[1].set_xlim(-10,10)
                    ax[2].set_ylabel('A'        )
                    #ax[2].set_xlim(-10,10)
                    #ax[2].set_ylim(0,0.4)
                    ax[3].set_ylabel('ImG'        )
                    ax[3].set_xlabel('Frequency')
                    #ax[3].set_ylim(-0.4,0.4)
                    #ax[3].set_xlim(-10,10)
                    #
                    for iax in range(0,4):
                        if bool(xylabel[iax][0]) :
                            ax[iax].set_xlabel = xylabel[iax][0]
                        if bool(xylabel[iax][1]) :
                            ax[iax].set_ylabel = xylabel[iax][1]
                        if bool(xylim[iax][0]) :
                            ax[iax].xlim = xylim[iax][0]
                        if bool(xylim[iax][1]) :
                            ax[iax].ylim = xylim[iax][1]
                    #
                    for delta in self.listdelta:
                        #
                        for nx in range(self.xstart,self.xend,self.xstep):
                            #
                            file_arr_s = np.loadtxt(self.dirfile+'/'+
                                                 classfile[0]+'_vs_w_k_'+str(k)+'_iter_'+str(iiter)+'_delta_'+str(delta)+'.dat')
                            file_arr_A = np.loadtxt(self.dirfile+'/'+
                                                 classfile[1]+'_vs_w_k_'+str(k)+'_iter_'+str(iiter)+'_delta_'+str(delta)+'.dat')
                            #
                            x = file_arr_s[:,1]
                            Re_S = file_arr_s[:,2]
                            Re_S_fit = file_arr_s[:,4]
                            Im_S = file_arr_s[:,3]
                            Im_S_fit = file_arr_s[:,5]
                            A    = file_arr_A[:,4]
                            A_fit= file_arr_A[:,7]
                            G= file_arr_A[:,3]
                            G_fit= file_arr_A[:,6]
                            #
                            ax[0].set_title("Self-energy at k="+str(k)+" for rs=4. qpoints = "+str(nq),y=1.05)
                            #
                            xstr = 'nx=' + str(nx) 
                            #
                            ax[0].plot(x,Re_S,linestyle='-',linewidth=0.5,label=xstr)
                            ax[1].plot(x,Im_S,linestyle='-',linewidth=0.5,label=xstr)
                            ax[0].plot(x,Re_S_fit,linestyle='-',linewidth=0.5,label=xstr)
                            ax[1].plot(x,Im_S_fit,linestyle='-',linewidth=0.5,label=xstr)
                            ax[2].plot(x,A,marker='.',label=xstr)
                            ax[2].plot(x,A_fit,linestyle='-',linewidth=0.5,label='fitted')
                            ax[3].plot(x,G,marker='.',label=xstr)
                            ax[3].plot(x,G_fit,linestyle='-',linewidth=0.5,label='fitted')
                            leg1 = ax[0].legend(loc=4,prop={'size': 6}) 
                            leg2 = ax[1].legend(loc=4,prop={'size': 6})
                            leg3 = ax[2].legend(loc=4,prop={'size': 6})
                            leg4 = ax[3].legend(loc=4,prop={'size': 6})
                            #
                            plt.subplots_adjust(hspace=0.3)
                            #
                            #
                        plt.savefig(fileout+'_k_'+str(k)+'_nq_'+str(nq)+'_iter_'+str(iiter)+'_delta_'+str(delta)+'.png',dpi=300)
                        plt.close()
            

    def pol_and_wscreen(self,xylabel,xylim,fileout='conv_tests',classfile=['polarizability_plot','wscreen_plot']):
        """
        this functions plots the polarizability and the W screened potential 
        """
        for q in self.listpps:
            #
            for iiter in self.listiter:
                #
                for nk in range(self.vstart,self.vend,self.vstep):
                    #
                    fig = plt.figure()
                    #
                    ax = []
                    for iax in range(0,4):
                        ax.append(fig.add_subplot(4,1,iax+1))
                    #
                    ax[0].set_ylabel('Re(P)')
                    ax[0].set_xlim(-10,10)
                    ax[1].set_ylabel('Im(P)')
                    ax[1].set_xlim(-10,10)
                    ax[2].set_ylabel('Re(Wcorr)' )
                    ax[2].set_xlim(-10,10)
                    ax[3].set_ylabel('Im(Wcorr)' )
                    ax[3].set_xlabel('Frequency')
                    ax[3].set_xlim(-10,10)
                    #
                    for iax in range(0,4):
                        if bool(xylabel[iax][0]) :
                            ax[iax].set_xlabel = xylabel[iax][0]
                        if bool(xylabel[iax][1]) :
                            ax[iax].set_ylabel = xylabel[iax][1]
                        if bool(xylim[iax][0]) :
                            ax[iax].xlim = xylim[iax][0]
                        if bool(xylim[iax][1]) :
                            ax[iax].ylim = xylim[iax][1]
                    #
                    for delta in self.listdelta:
                        #
                        for nx in range(self.xstart,self.xend,self.xstep):
                            #
                            file_arr_P = np.loadtxt(self.dirfile+'/'+
                                                 classfile[0]+'_vs_w_q_'+str(q)+'_iter_'+str(iiter)+'_delta_'+str(delta)+'.dat')
                            file_arr_W = np.loadtxt(self.dirfile+'/'+
                                                 classfile[1]+'_vs_w_q_'+str(q)+'_iter_'+str(iiter)+'_delta_'+str(delta)+'.dat')
                            #
                            x = file_arr_P[:,0]
                            Re_P = file_arr_P[:,1]
                            Re_Pfit = file_arr_P[:,5]
                            Im_P = file_arr_P[:,2]
                            Im_Pfit = file_arr_P[:,6]
                            Re_Wcorr = file_arr_W[:,1]
                            Re_Wcorrfit = file_arr_W[:,5]
                            Im_Wcorr = file_arr_W[:,2]
                            Im_Wcorrfit = file_arr_W[:,6]
                            #
                            #
                            ax[0].set_title("Polarizability at q="+str(q)+" for rs=4. qpoints = "+str(nk),y=1.05)
                            #
                            xstr = 'nxh=' + str(nx) 
                            #
                            ax[0].plot(x,Re_P,marker='.',label=xstr)
                            ax[0].plot(x,Re_Pfit,linestyle='-',linewidth=0.5,label='fitted')
                            ax[1].plot(x,Im_P,marker='.',label=xstr)
                            ax[1].plot(x,Im_Pfit,linestyle='-',linewidth=0.5,label='fitted')
                            ax[2].plot(x,Re_Wcorr,marker='.',label=xstr)
                            ax[2].plot(x,Re_Wcorrfit,linestyle='-',linewidth=0.5,label='fitted')
                            ax[3].plot(x,Im_Wcorr,marker='.',label=xstr)
                            ax[3].plot(x,Im_Wcorrfit,linestyle='-',linewidth=0.5,label='fitted')
                            leg1 = ax[0].legend(loc=4,prop={'size': 6}) 
                            leg2 = ax[1].legend(loc=4,prop={'size': 6})
                            leg3 = ax[2].legend(loc=4,prop={'size': 6})
                            leg4 = ax[3].legend(loc=4,prop={'size': 6})
                            #
                            plt.subplots_adjust(hspace=0.3)
                            #
                            #
                        plt.savefig(fileout+'_q_'+str(q)+'_nk_'+str(nk)+'_iter_'+str(iiter)+'_delta_'+str(delta)+'.png',dpi=300)
                        plt.close()
            

    def sc_sigma_and_greenf(self,xylabel,xylim,fileout='compar',classfile=['sigma_plot','greenfunc_plot'],W0=False):
        """
        this functions plots the self-energy and the greensfunction comparing the +1 iteration with the current one
        """
        for k in self.listpps:
            #
            for iiter in self.listiter:
                #
                for nq in range(self.vstart,self.vend,self.vstep):
                    #
                    fig = plt.figure()
                    #
                    ax = []
                    for iax in range(0,4):
                        ax.append(fig.add_subplot(4,1,iax+1))
                    #
                    ax[0].set_ylabel('Re(Sigma)')
                    # ax[0].set_xlim(-10,10)
                    ax[1].set_ylabel('Im(Sigma)')
                    # ax[1].set_xlim(-10,10)
                    ax[2].set_ylabel('A'        )
                    #ax[2].set_xlim(-10,10)
                    #ax[2].set_ylim(0,5)
                    ax[3].set_ylabel('A'        )
                    ax[3].set_xlabel('Frequency')
                    #ax[3].set_ylim(0,0.5)
                    #ax[3].set_xlim(-10,10)
                    #
                    for iax in range(0,4):
                        if bool(xylabel[iax][0]) :
                            ax[iax].set_xlabel = xylabel[iax][0]
                        if bool(xylabel[iax][1]) :
                            ax[iax].set_ylabel = xylabel[iax][1]
                        if bool(xylim[iax][0]) :
                            ax[iax].xlim = xylim[iax][0]
                        if bool(xylim[iax][1]) :
                            ax[iax].ylim = xylim[iax][1]
                    #
                    for delta in self.listdelta:
                        #
                        for nx in range(self.xstart,self.xend,self.xstep):
                            #
                            file_arr_s0 = np.loadtxt(self.dirfile+'/'+
                                                 classfile[0]+'_vs_w_k_'+str(k)+'_iter_'+str(iiter)+'_delta_'+str(delta)+'.dat')
                            file_arr_A0 = np.loadtxt(self.dirfile+'/'+
                                                 classfile[1]+'_vs_w_k_'+str(k)+'_iter_'+str(iiter)+'_delta_'+str(delta)+'.dat')
                            file_arr_s = np.loadtxt(self.dirfile+'/'+
                                                 classfile[0]+'_vs_w_k_'+str(k)+'_iter_'+str(iiter+1)+'_delta_'+str(delta)+'.dat')
                            file_arr_A = np.loadtxt(self.dirfile+'/'+
                                                 classfile[1]+'_vs_w_k_'+str(k)+'_iter_'+str(iiter+1)+'_delta_'+str(delta)+'.dat')
                            #
                            x = file_arr_s[:,1]
                            x0 = file_arr_s0[:,1]
                            Re_S0 = file_arr_s0[:,2]
                            Im_S0 = file_arr_s0[:,3]
                            A0 = file_arr_A0[:,4]
                            A0_fit = file_arr_A0[:,7]
                            Re_S = file_arr_s[:,2]
                            Im_S = file_arr_s[:,3]
                            A    = file_arr_A[:,4]
                            A_fit= file_arr_A[:,7]
                            #
                            if(W0):
                                ax[0].set_title("G"+str(iiter)+"W0 vs G"+str(iiter+1)+"W0 at k="+str(k)+" for rs=4. qpoints = "+str(nq),y=1.05)
                                xstr0 = 'G'+str(iiter)+'W0'
                                xstr1 = 'G'+str(iiter+1)+'W0'
                            else:
                                ax[0].set_title("GW"+str(iiter)+" vs GW"+str(iiter+1)+" at k="+str(k)+" for rs=4. qpoints = "+str(nq),y=1.05)
                                xstr0 = "GW"+str(iiter) 
                                xstr1 = "GW"+str(iiter+1) 
                            #
                            #
                            ax[0].plot(x0,Re_S0,linestyle='-',linewidth=0.5,label=xstr0)
                            ax[1].plot(x0,Im_S0,linestyle='-',linewidth=0.5,label=xstr0)
                            ax[2].plot(x0,A0_fit,linestyle='-',linewidth=0.5,label=xstr0)
                            ax[3].plot(x0,A0_fit,linestyle='-',linewidth=0.5,label=xstr0)
                            #
                            #
                            ax[0].plot(x,Re_S,linestyle='-',linewidth=0.5,label=xstr1)
                            ax[1].plot(x,Im_S,linestyle='-',linewidth=0.5,label=xstr1)
                            # ax[2].plot(x,A,marker='.',label=xstr1)
                            ax[2].plot(x,A_fit,linestyle='-',linewidth=0.5,label=xstr1)
                            ax[3].plot(x,A_fit,linestyle='-',linewidth=0.5,label=xstr1)
                            leg1 = ax[0].legend(loc=4,prop={'size': 6}) 
                            leg2 = ax[1].legend(loc=4,prop={'size': 6})
                            leg3 = ax[2].legend(loc=4,prop={'size': 6})
                            leg4 = ax[3].legend(loc=4,prop={'size': 6})
                            #
                            plt.subplots_adjust(hspace=0.3)
                            #
                            #
                        plt.savefig(fileout+'_k_'+str(k)+'_nq_'+str(nq)+'_iter_'+str(iiter)+'_delta_'+str(delta)+'.png',dpi=300)
                        plt.close()


    def sc_pol_and_wscreen(self,xylabel,xylim,fileout='conv_tests',classfile=['polarizability_plot','wscreen_plot']):
        """
        this functions plots the polarizability and the W screened potential comparing the +1 iteration with the current one 
        """
        for q in self.listpps:
            #
            for iiter in self.listiter:
                #
                for nk in range(self.vstart,self.vend,self.vstep):
                    #
                    fig = plt.figure()
                    #
                    ax = []
                    for iax in range(0,4):
                        ax.append(fig.add_subplot(4,1,iax+1))
                    #
                    ax[0].set_ylabel('Re(P)')
                    ax[0].set_xlim(-10,10)
                    ax[1].set_ylabel('Im(P)')
                    ax[1].set_xlim(-10,10)
                    ax[2].set_ylabel('Re(Wcorr)' )
                    ax[2].set_xlim(-10,10)
                    ax[3].set_ylabel('Im(Wcorr)' )
                    ax[3].set_xlabel('Frequency')
                    ax[3].set_xlim(-10,10)
                    #
                    for iax in range(0,4):
                        if bool(xylabel[iax][0]) :
                            ax[iax].set_xlabel = xylabel[iax][0]
                        if bool(xylabel[iax][1]) :
                            ax[iax].set_ylabel = xylabel[iax][1]
                        if bool(xylim[iax][0]) :
                            ax[iax].xlim = xylim[iax][0]
                        if bool(xylim[iax][1]) :
                            ax[iax].ylim = xylim[iax][1]
                    #
                    for delta in self.listdelta:
                        #
                        for nx in range(self.xstart,self.xend,self.xstep):
                            #
                            file_arr_P0 = np.loadtxt(self.dirfile+'/'+
                                                  classfile[0]+'_vs_w_q_'+str(q)+'_iter_'+str(iiter)+'_delta_'+str(delta)+'.dat')
                            file_arr_W0 = np.loadtxt(self.dirfile+'/'+
                                                  classfile[1]+'_vs_w_q_'+str(q)+'_iter_'+str(iiter)+'_delta_'+str(delta)+'.dat')
                            #
                            file_arr_P = np.loadtxt(self.dirfile+'/'+
                                                  classfile[0]+'_vs_w_q_'+str(q)+'_iter_'+str(iiter+1)+'_delta_'+str(delta)+'.dat')
                            file_arr_W = np.loadtxt(self.dirfile+'/'+
                                                  classfile[1]+'_vs_w_q_'+str(q)+'_iter_'+str(iiter+1)+'_delta_'+str(delta)+'.dat')
                            #
                            x = file_arr_P[:,0]
                            Re_P = file_arr_P[:,1]
                            Re_Pfit = file_arr_P[:,5]
                            Im_P = file_arr_P[:,2]
                            Im_Pfit = file_arr_P[:,6]
                            Re_Wcorr = file_arr_W[:,1]
                            Re_Wcorrfit = file_arr_W[:,5]
                            Im_Wcorr = file_arr_W[:,2]
                            Im_Wcorrfit = file_arr_W[:,6]
                            #
                            x0 = file_arr_P0[:,0]
                            Re_P0 = file_arr_P0[:,1]
                            Re_Pfit0 = file_arr_P0[:,5]
                            Im_P0 = file_arr_P0[:,2]
                            Im_Pfit0 = file_arr_P0[:,6]
                            Re_Wcorr0 = file_arr_W0[:,1]
                            Re_Wcorrfit0 = file_arr_W0[:,5]
                            Im_Wcorr0 = file_arr_W0[:,2]
                            Im_Wcorrfit0 = file_arr_W0[:,6]
                            #
                            ax[0].set_title("Polarizability at iter "+str(iiter)+" vs iter "+str(iiter+1)+" at q="+str(q)+" for rs=4. kpoints = "+str(nk),y=1.05)
                            #
                            xstr = 'iter '+str(iiter)
                            #
                            ax[0].plot(x0,Re_P0,marker='.',label=xstr)
                            ax[0].plot(x0,Re_Pfit0,linestyle='-',linewidth=0.5,label='fitted')
                            ax[1].plot(x0,Im_P0,marker='.',label=xstr)
                            ax[1].plot(x0,Im_Pfit0,linestyle='-',linewidth=0.5,label='fitted')
                            ax[2].plot(x0,Re_Wcorr0,marker='.',label=xstr)
                            ax[2].plot(x0,Re_Wcorrfit0,linestyle='-',linewidth=0.5,label='fitted')
                            ax[3].plot(x0,Im_Wcorr0,marker='.',label=xstr)
                            ax[3].plot(x0,Im_Wcorrfit0,linestyle='-',linewidth=0.5,label='fitted')
                            #
                            xstr = 'iter '+str(iiter+1)
                            #
                            ax[0].plot(x,Re_P,marker='.',label=xstr)
                            ax[0].plot(x,Re_P,linestyle='-',linewidth=0.5,label='fitted')
                            ax[0].plot(x,Re_Pfit,linestyle='-',linewidth=0.5,label='fitted')
                            ax[1].plot(x,Im_P,marker='.',label=xstr)
                            ax[1].plot(x,Im_Pfit,linestyle='-',linewidth=0.5,label='fitted')
                            ax[2].plot(x,Re_Wcorr,marker='.',label=xstr)
                            ax[2].plot(x,Re_Wcorrfit,linestyle='-',linewidth=0.5,label='fitted')
                            ax[3].plot(x,Im_Wcorr,marker='.',label=xstr)
                            ax[3].plot(x,Im_Wcorrfit,linestyle='-',linewidth=0.5,label='fitted')
                            #
                            leg1 = ax[0].legend(loc=4,prop={'size': 6}) 
                            leg2 = ax[1].legend(loc=4,prop={'size': 6})
                            leg3 = ax[2].legend(loc=4,prop={'size': 6})
                            leg4 = ax[3].legend(loc=4,prop={'size': 6})
                            #
                            plt.subplots_adjust(hspace=0.3)
                            #
                            #
                        plt.savefig(fileout+'_q_'+str(q)+'_nk_'+str(nk)+'_iter_'+str(iiter)+'_delta_'+str(delta)+'.png',dpi=300)
                        plt.close()
            

    def sc_occ_number(self,xylabel,xylim,fileout='sc_occ_number',classfile=['occ_number']):
        """
        this functions plots the sc occupation number comparing the +1 iteration with the current one
        """
        for iiter in self.listiter:
            #
            for nq in range(self.vstart,self.vend,self.vstep):
                #
                fig = plt.figure()
                #
                ax = []
                for iax in range(0,1):
                    ax.append(fig.add_subplot(1,1,iax+1))
                #
                ax[0] = fig.add_subplot(111)
                ax[0].set_ylabel('n_k')
                ax[0].set_ylabel('k')
                ax[0].set_xlim(-0.1,4)
                #
                for iax in range(0,4):
                    if bool(xylabel[iax][0]) :
                        ax[iax].set_xlabel = xylabel[iax][0]
                    if bool(xylabel[iax][1]) :
                        ax[iax].set_ylabel = xylabel[iax][1]
                    if bool(xylim[iax][0]) :
                        ax[iax].xlim = xylim[iax][0]
                    if bool(xylim[iax][1]) :
                        ax[iax].ylim = xylim[iax][1]
                #
                for delta in self.listdelta:
                    #
                    for nx in range(self.xstart,self.xend,self.xstep):
                        #
                        file_arr_s = np.loadtxt(self.dirfile+'/'+
                                             str(iiter+1)+'_'+classfile[0]+'_delta_'+str(delta)+'.dat')
                        file_arr_s0 = np.loadtxt(self.dirfile+'/'+
                             str(iiter)+'_'+classfile[0]+'_delta_'+str(delta)+'.dat')
                        #
                        k = file_arr_s[:,0]
                        k0 = file_arr_s0[:,0]
                        occ_k = file_arr_s[:,1]
                        occ_k0 = file_arr_s0[:,1]
                        #
                        ax[0].set_title("Occupation number in GW0. "+str(iiter)+"_n_k vs "+str(iiter+1)+"_n_k. \n All units are in [a_f].")
                        #
                        labeln = str(iiter+1)+"_n_k" 
                        #
                        ax[0].plot(k,occ_k,linestyle='-',linewidth=0.5,label=labeln)
                        #
                        label0 = str(iiter)+"_n_k" 
                        #
                        ax[0].plot(k0,occ_k0,linestyle='-',linewidth=0.5,label=label0)
                        #
                        #
                    plt.savefig(fileout+'_nq_'+str(nq)+'_iter_'+str(iiter)+'_delta_'+str(delta)+'.png',dpi=300)
                    plt.close()
        

    def occ_number(self,xylabel,xylim,fileout='occ_number',classfile=['occ_number']):
        """
        this functions plots the occupation number
        """
        for iiter in self.listiter:
            #
            for nq in range(self.vstart,self.vend,self.vstep):
                #
                fig = plt.figure()
                #
                ax = []
                for iax in range(0,1):
                    ax.append(fig.add_subplot(1,1,iax+1))
                #
                ax[0].set_ylabel('n_k')
                ax[0].set_xlabel('k')
                ax[0].set_xlim(-0.1,4)
                #
                for iax in range(0,4):
                    if bool(xylabel[iax][0]) :
                        ax[iax].set_xlabel = xylabel[iax][0]
                    if bool(xylabel[iax][1]) :
                        ax[iax].set_ylabel = xylabel[iax][1]
                    if bool(xylim[iax][0]) :
                        ax[iax].xlim = xylim[iax][0]
                    if bool(xylim[iax][1]) :
                        ax[iax].ylim = xylim[iax][1]
                #
                for delta in self.listdelta:
                    #
                    for nx in range(self.xstart,self.xend,self.xstep):
                        #
                        file_arr_s = np.loadtxt(self.dirfile+'/'+
                                             str(iiter)+'_'+classfile[0]+'_delta_'+str(delta)+'.dat')
                        #
                        k = file_arr_s[:,0]
                        occ_k = file_arr_s[:,1]
                        #
                        ax[0].set_title("Iteration n. "+str(iiter)+". Occupation number. \n All units are in [a_f].")
                        #
                        ax[0].plot(k,occ_k,linestyle='-',linewidth=0.5,label='n_k')
                        leg1 = ax[0].legend(loc=4,prop={'size': 6}) 
                        #
                        #
                    plt.savefig(fileout+'_nq_'+str(nq)+'_iter_'+str(iiter)+'_delta_'+str(delta)+'.png',dpi=300)
                    plt.close()
        
    def sc_dos(self,xylabel,xylim,fileout='sc_dos',classfile=['spectrfunc']):
        """
        this functions plots the sc occupation number comparing the +1 iteration with the current one
        """
        for iiter in self.listiter:
            #
            for nq in range(self.vstart,self.vend,self.vstep):
                #
                fig = plt.figure()
                #
                ax = []
                for iax in range(0,1):
                    ax.append(fig.add_subplot(1,1,iax+1))
                #
                ax[0].set_ylabel('DOS')
                ax[0].set_xlabel('omega')
                ax[0].set_xlim(-5,5)
                ax[0].set_ylim(0,0.13)
                #
                for iax in range(0,4):
                    if bool(xylabel[iax][0]) :
                        ax[iax].set_xlabel = xylabel[iax][0]
                    if bool(xylabel[iax][1]) :
                        ax[iax].set_ylabel = xylabel[iax][1]
                    if bool(xylim[iax][0]) :
                        ax[iax].xlim = xylim[iax][0]
                    if bool(xylim[iax][1]) :
                        ax[iax].ylim = xylim[iax][1]
                #
                for delta in self.listdelta:
                    #
                    for nx in range(self.xstart,self.xend,self.xstep):
                        #
                        file_arr_s = np.loadtxt(self.dirfile+'/'+
                                             str(iiter+1)+'_'+classfile[0]+'_delta_'+str(delta)+'.dat')
                        file_arr_s0 = np.loadtxt(self.dirfile+'/'+
                             str(iiter)+'_'+classfile[0]+'_delta_'+str(delta)+'.dat')
                        #
                        omega = file_arr_s[:,0]
                        omega0 = file_arr_s0[:,0]
                        f_omega = file_arr_s[:,1]
                        f_omega0 = file_arr_s0[:,1]
                        #
                        ax[0].set_title("Spectral function in GW"+str(iiter)+" vs GW"+str(iiter+1)+". \n All units are in [fermi].")
                        #
                        labeln = "iter "+str(iiter+1) 
                        #
                        ax[0].plot(omega,f_omega,linestyle='-',linewidth=0.5,label=labeln)
                        #
                        label0 = "iter "+str(iiter) 
                        #
                        ax[0].plot(omega0,f_omega0,linestyle='-',linewidth=0.5,label=label0)
                        #
                        leg1 = ax[0].legend(loc=4,prop={'size': 6}) 
                        #
                        #
                    plt.savefig(fileout+'_nq_'+str(nq)+'_iter_'+str(iiter)+'_delta_'+str(delta)+'.png',dpi=300)
                    plt.close()
        

    def dos(self,xylabel,xylim,fileout='dos',classfile=['spectrfunc']):
        """
        this functions plots the occupation number
        """
        for iiter in self.listiter:
            #
            for nq in range(self.vstart,self.vend,self.vstep):
                #
                fig = plt.figure()
                #
                ax = []
                for iax in range(0,1):
                    ax.append(fig.add_subplot(1,1,iax+1))
                #
                ax[0].set_ylabel('DOS')
                ax[0].set_xlabel('omega')
                ax[0].set_xlim(-5,5)
                ax[0].set_ylim(0,0.13)
                #
                for iax in range(0,4):
                    if bool(xylabel[iax][0]) :
                        ax[iax].set_xlabel = xylabel[iax][0]
                    if bool(xylabel[iax][1]) :
                        ax[iax].set_ylabel = xylabel[iax][1]
                    if bool(xylim[iax][0]) :
                        ax[iax].xlim = xylim[iax][0]
                    if bool(xylim[iax][1]) :
                        ax[iax].ylim = xylim[iax][1]
                #
                for delta in self.listdelta:
                    #
                    for nx in range(self.xstart,self.xend,self.xstep):
                        #
                        file_arr_s = np.loadtxt(self.dirfile+'/'+
                                             str(iiter)+'_'+classfile[0]+'_delta_'+str(delta)+'.dat')
                        #
                        omega = file_arr_s[:,0]
                        f_omega = file_arr_s[:,1]
                        #
                        ax[0].set_title("Iteration n. "+str(iiter)+". Spectral function. \n All units are in [fermi].")
                        #
                        ax[0].plot(omega,f_omega,linestyle='-',linewidth=0.5,label='DOS')
                        leg1 = ax[0].legend(loc=4,prop={'size': 6}) 
                        #
                        #
                    plt.savefig(fileout+'_nq_'+str(nq)+'_iter_'+str(iiter)+'_delta_'+str(delta)+'.png',dpi=300)
                    plt.close()
        

    def superimpose_dos(self,xylabel,xylim,fileout='dos',classfile=['spectrfunc']):
        """
        this functions plots the occupation number
        """
        #
        #
        fig = plt.figure()
        #
        ax = []
        for iax in range(0,1):
            ax.append(fig.add_subplot(1,1,iax+1))
        #
        for delta in self.listdelta:
            #
            ax[0].set_ylabel('DOS')
            ax[0].set_xlabel('omega')
            ax[0].set_xlim(-1,1)
            ax[0].set_ylim(0,1.5)
            #
            ax[0].set_title("Spectral function in SC-GW calculations. \n All units are in [Ha].")
            #
            file_arr_s = np.loadtxt(self.dirfile+'/'+
                                 'free'+'_'+classfile[0]+'_delta_'+str(delta)+'.dat')
            #
            omega = file_arr_s[:,0]
            f_omega = file_arr_s[:,1]
            labell='free-theoretical'
            #
            ax[0].plot(omega,f_omega,linestyle='-',linewidth=2.0,label=labell)
            #
            for iiter in self.listiter:
                #
                file_arr_s = np.loadtxt(self.dirfile+'/'+
                                     str(iiter)+'_'+classfile[0]+'_delta_'+str(delta)+'.dat')
                #
                omega = file_arr_s[:,0]
                f_omega = file_arr_s[:,1]
                if(iiter==-1):
                    labell='free'
                elif(iiter==0):
                    labell='HF'
                elif(iiter==1):
                    labell='G1W1'
                elif(iiter==8):
                    labell='sc-GW'
                #
                ax[0].plot(omega,f_omega,linestyle='-',linewidth=2.0,label=labell)
                #
            leg1 = ax[0].legend(loc=4,prop={'size': 6}) 
            plt.savefig(fileout+'_superimpose_delta_'+str(delta)+'.png',dpi=300)
            plt.close()
        
    def superimpose_spectralpot(self,xylabel,xylim,fileout='spectrpot',classfile=['spectrpot']):
        """
        this functions plots the occupation number
        """
        #
        #
        fig = plt.figure()
        #
        ax = []
        for iax in range(0,1):
            ax.append(fig.add_subplot(1,1,iax+1))
        #
        for delta in self.listdelta:
            #
            ax[0].set_ylabel('V')
            ax[0].set_xlabel('omega')
            ax[0].set_xlim(-1,1)
            ax[0].set_ylim(-0.35,0.35)
            #
            ax[0].set_title("Preliminary results for spectral potential of [6] in SC-GW calculations. \n All units are in [Ha].")
            #
            file_arr_s = np.loadtxt(self.dirfile+'/'+
                                 'free'+'_'+classfile[0]+'_delta_'+str(delta)+'.dat')
            #
            omega = file_arr_s[:,0]
            f_omega = file_arr_s[:,1]
            labell='free theoretical'
            #
            ax[0].plot(omega,f_omega,linestyle='-',linewidth=2.0,label=labell)
            for iiter in self.listiter:
                #
                file_arr_s = np.loadtxt(self.dirfile+'/'+
                                     str(iiter)+'_'+classfile[0]+'_delta_'+str(delta)+'.dat')
                #
                omega = file_arr_s[:,0]
                f_omega = file_arr_s[:,1]
                if(iiter==-1):
                    labell='free'
                elif(iiter==0):
                    labell='HF'
                elif(iiter==1):
                    labell='G1W1'
                elif(iiter==8):
                    labell='sc-GW'
                #
                ax[0].plot(omega,f_omega,linestyle='-',linewidth=2.0,label=labell)
                #
            leg1 = ax[0].legend(loc=4,prop={'size': 6}) 
            plt.savefig(fileout+'_superimpose_delta_'+str(delta)+'.png',dpi=300)
            plt.close()
        
    def superimpose_sigma_greenf(self,xylabel,xylim,fileout='sigma_greenf_',classfile=['sigma_plot','greenfunc_plot'],W0=False):
        """
        this functions plots the self-energy and the greensfunction comparing the +1 iteration with the current one
        """
        for k in self.listpps:
            #
            fig = plt.figure()
            ax = []
            ax.append(fig.add_subplot(4,1,1))
            ax.append(fig.add_subplot(4,1,2))
            ax.append(fig.add_subplot(4,1,3))
            ax.append(fig.add_subplot(4,1,4))
            #
            ax[0].set_ylabel('Re(Sigma)')
            ax[1].set_ylabel('Im(Sigma)')
            ax[2].set_ylabel('A'        )
            ax[0].set_xlim(-10,10)
            ax[1].set_xlim(-10,10)
            ax[2].set_xlim(-10,10)
            ax[2].set_ylim(0,5)
            ax[3].set_ylabel('A'        )
            ax[3].set_xlabel('Frequency')
            ax[3].set_ylim(0,0.5)
            ax[3].set_xlim(-10,10)
            ax[0].set_title("G"+str(self.listiter[0]-1)+"W"+str(self.listiter[0]-1)+" vs G"+str(self.listiter[1]-1)+"W"+str(self.listiter[1]-1)+" at k="+str(k)+" for rs=4.",y=1.05)
            #
            ax[0].tick_params(
            axis='x',          # changes apply to the x-axis
            which='both',      # both major and minor ticks are affected
            bottom=False,      # ticks along the bottom edge are off
            top=False,         # ticks along the top edge are off
            labelbottom=False) # labels along the bottom edge are off
            #
            ax[1].tick_params(
            axis='x',          # changes apply to the x-axis
            which='both',      # both major and minor ticks are affected
            bottom=False,      # ticks along the bottom edge are off
            top=False,         # ticks along the top edge are off
            labelbottom=False) # labels along the bottom edge are off
            #
            ax[2].tick_params(
            axis='x',          # changes apply to the x-axis
            which='both',      # both major and minor ticks are affected
            bottom=False,      # ticks along the bottom edge are off
            top=False,         # ticks along the top edge are off
            labelbottom=False) # labels along the bottom edge are off
            #
            for iiter in self.listiter:
                #
                for delta in self.listdelta:
                    #
                    for nx in range(self.xstart,self.xend,self.xstep):
                        #
                        file_arr_s0 = np.loadtxt(self.dirfile+'/'+
                                             classfile[0]+'_vs_w_k_'+str(k)+'_iter_'+str(iiter)+'_delta_'+str(delta)+'.dat')
                        file_arr_A0 = np.loadtxt(self.dirfile+'/'+
                                             classfile[1]+'_vs_w_k_'+str(k)+'_iter_'+str(iiter)+'_delta_'+str(delta)+'.dat')
                        #
                        x0 = file_arr_s0[:,0]
                        Re_S0 = file_arr_s0[:,1]
                        Im_S0 = file_arr_s0[:,2]
                        A0 = file_arr_A0[:,3]
                        A0_fit = file_arr_A0[:,6]
                        #
                        if(iiter==-1):
                            labell='free'
                        elif(iiter==0):
                            labell='HF'
                        elif(iiter==1):
                            labell='G0W0'
                        elif(iiter==8):
                            labell='sc-GW'
                        #
                        #
                        ax[0].plot(x0,Re_S0,linestyle='-',linewidth=2.0,label=labell)
                        ax[1].plot(x0,Im_S0,linestyle='-',linewidth=2.0,label=labell)
                        ax[2].plot(x0,A0_fit,linestyle='-',linewidth=2.0,label=labell)
                        ax[3].plot(x0,A0_fit,linestyle='-',linewidth=2.0,label=labell)
                        #
                        #
                        leg1 = ax[0].legend(loc=4,prop={'size': 6}) 
                        leg2 = ax[1].legend(loc=4,prop={'size': 6})
                        leg3 = ax[2].legend(loc=4,prop={'size': 6})
                        leg4 = ax[3].legend(loc=4,prop={'size': 6})
                        #
                        plt.subplots_adjust(hspace=0)
                        #
                        #
            plt.savefig(fileout+'superimpose'+'_k_'+str(k)+'_delta_'+str(delta)+'.png',dpi=300)
            plt.show()
            plt.close()
