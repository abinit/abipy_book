import numpy as np
import math 
import matplotlib.pyplot as plt
from math import factorial
from matplotlib.pyplot import cm
from scipy import special

from ipywidgets import interact, interactive, fixed, interact_manual
import ipywidgets as widgets
import plotly.io as pio
pio.templates.default = "none"
from plotly.subplots import make_subplots
import plotly.graph_objects as go



def psi(Q,n,alpha):
    
    N = np.sqrt(1/(2**n*factorial(n))) * (alpha/np.pi)**(1/4)
    psi= N * special.eval_hermite(n,np.sqrt(alpha)*Q) * np.exp(-(alpha*Q**2)/2)
    
    return psi

def energy(n,omega):
    E=(n+0.5)*omega
    return E

def potential(omega,Q):
    return (1/2)*omega**2*Q**2
    
def inverse_potential(E,omega):
    return np.sqrt((2*E))/omega

def get_S(omega,Q_offset):
    return omega*Q_offset**2/2

def get_omega(mu,k):
    return np.sqrt(k/mu)

def lineshape(n,Energy,omega,Q_offset,E_zpl,broadening=0.1):
    sigma=broadening*omega
    A=np.zeros(len(Energy))
    S=get_S(omega,Q_offset)
    for i in range(n):
        A+= np.exp(-S)*(S**i)/(factorial(i)) * (1/sigma*np.sqrt(2*np.pi))*np.exp(-0.5*(((Energy-(E_zpl-omega*i))/sigma)**2))

    return A/max(A)

def lineshape_envelope(n,Energy,omega,Q_offset,E_zpl,broadening=0.7):
    
    sigma_envelop=omega*broadening
    A_envelop=np.zeros(len(Energy))
    S=get_S(omega,Q_offset)        
    for i in range(n):
        A_envelop+=np.exp(-S)*(S**i)/(factorial(i)) * (1/sigma_envelop*np.sqrt(2*np.pi))*np.exp(-0.5*(((Energy-(E_zpl-omega*i))/sigma_envelop)**2))
    
    return A_envelop/max(A_envelop)

def get_plotly_lineshape(slider="frequency",num_step=30,range_k=[0.2,2],range_Q=[0.5,3],default_k=0.8,default_Q=2,mu=0.5):
    """
    slider= "frequency" or "Delta_Q"
    """
    ## parameters, might be changed
    mu=mu
    E_zpl=15
    Energy=np.linspace(E_zpl-30,E_zpl+2,10000)
    width_paras=4
    n_vib=5
    
    if slider == "frequency":
        ks=np.linspace(range_k[0],range_k[1],num_step)
        Q_offsets=np.array([default_Q]*len(ks))
    if slider == "Delta_Q":
        Q_offsets=np.linspace(range_Q[0],range_Q[1],num_step)
        ks=np.array([default_k]*len(Q_offsets))

    ## do not touch
    
    omegas=[get_omega(mu,k) for k in ks]
    
    As=[lineshape(n=20,Energy=Energy,omega=omega,Q_offset=Q_offset,E_zpl=E_zpl) 
        for omega,Q_offset in zip(omegas,Q_offsets)]
    
    As_envelope=[lineshape_envelope(n=20,Energy=Energy,omega=omega,Q_offset=Q_offset,E_zpl=E_zpl) 
                 for omega,Q_offset in zip(omegas,Q_offsets)]
    
    Q_para_1=[np.linspace(-width_paras,width_paras,1000) for omega,Q_offset in zip(omegas,Q_offsets)]
    
    pot_1s=[potential(omega,Q_para_1[i]) for i,omega in enumerate(omegas)]
    
    Q_para_2=[np.linspace(Q_offset-width_paras,Q_offset+width_paras,1000) for omega,Q_offset in zip(omegas,Q_offsets)]
    
#    Q_para_2=[np.linspace(Q_offset-width_paras,Q_offset+width_paras,1000) for Q_offset in Q_offsets]
    pot_2s=[(potential(omega,Q_para_2[i]-Q_offset)+E_zpl) for i,(omega,Q_offset) in enumerate(zip(omegas,Q_offsets))]
    

    ### Figure, do not touch
    
    fig = make_subplots(1, 2,
                            shared_yaxes=True,
                            horizontal_spacing=0.00,
                            column_widths=[0.7, 0.3],
                       )
    
    
    annotations_dict=[]
    
    for i in range(num_step):
        # lineshapes
        fig.add_trace(
            go.Scatter(
                visible=False,
                line=dict(color="black"),
                x=As[i],
                y=-Energy+E_zpl+0.5*omegas[i]),row=1, col=2)
    
        # lineshapes env
        fig.add_trace(
            go.Scatter(
                visible=False,
                line=dict(color="black"),
                x=As_envelope[i],
                y=-Energy+E_zpl+0.5*omegas[i]),row=1, col=2)
    
        # parabola gs
        fig.add_trace(
            go.Scatter(
                visible=False,
                line=dict(color="black"),
                x=Q_para_1[i],
                y=pot_1s[i]),row=1, col=1)
    
        # parabola_ex
        fig.add_trace(
            go.Scatter(
                visible=False,
                line=dict(color="black"),
                x=Q_para_2[i],
                y=pot_2s[i]),row=1, col=1)
    
    
        # anotations
        record_1 = go.layout.Annotation(
                text=f"&#969; = (k/&mu;)<sup>(1/2)</sup> = {np.round(omegas[i],2)}",
                y=20,
                x=0.5,
                showarrow=False,xref ='x2',yref ='y2',
            ) 
        record_2 = go.layout.Annotation(
                text=f"&#916;Q = {np.round(Q_offsets[i],2)}",
                y=22,
                x=0.5,
                showarrow=False,xref ='x2',yref ='y2',
            ) 
    
        record_3 = go.layout.Annotation(
                text=f"&mu; = {mu}",
                y=24,
                x=0.5,
                showarrow=False,xref ='x2',yref ='y2',
            ) 

        record_4 = go.layout.Annotation(
                text=f"S = {np.round(get_S(omegas[i],Q_offsets[i]),2)}",
                y=18,
                x=0.5,
                showarrow=False,xref ='x2',yref ='y2',
            ) 
        
        record_5 = go.layout.Annotation(
                text="(1/2)&#969;<sup>2</sup>Q<sup>2</sup>",
                y=-2,
                x=-2,
                showarrow=False,
            ) 
    
        record_6 = go.layout.Annotation(
                text="(1/2)&#969;<sup>2</sup>(Q-&#916;Q)<sup>2</sup>+E<sub>ZPL</sub>",
                y=E_zpl,
                x=Q_offsets[i]+4,
                showarrow=False,
            ) 
        annotations=[record_1,record_2,record_3,record_4,record_5,record_6]
        annotations_dict.append(annotations)
    
        # vib levels gs
        for j in range(n_vib):
            xmin=-inverse_potential(energy(j,omegas[i]),omegas[i])
            xmax=10
            fig.add_trace(
                go.Scatter(
                    visible=False,
                    x=[xmin,xmax],
                    y=[energy(j,omegas[i]),energy(j,omegas[i])],
                    line=dict(dash='dot',color="grey"),
                    mode="lines"),
                row=1, col=1)
    
        # vib levels ex
        for j in range(n_vib):
            xmin=-inverse_potential(energy(j,omegas[i]),omegas[i])+Q_offsets[i]
            xmax=+inverse_potential(energy(j,omegas[i]),omegas[i])+Q_offsets[i]
            fig.add_trace(
                go.Scatter(
                    visible=False,
                    x=[xmin,xmax],
                    y=[energy(j,omegas[i])+E_zpl,energy(j,omegas[i])+E_zpl],
                    line=dict(dash='dot',color="grey"),
                    mode="lines"),
                row=1, col=1)
    
    num_traces=4+2*n_vib
    
    # first image
    for i in range(0,num_traces):
        fig.data[i].visible=True
    fig.update_layout(
        annotations=annotations_dict[0],
    )

    # Create and add slider
    steps = []
    for index,i in enumerate(range(0, len(fig.data), num_traces)):
        if slider == "frequency":
            label_step=np.round(ks[index],2)
        if slider == "Delta_Q":
            label_step=np.round(Q_offsets[index],2)

        step = dict(
            method="update",
            args=[
                {"visible": [False] * len(fig.data)}, # update for traces
                {"title": "",    
                "annotations": annotations_dict[index]},],
            label=label_step,
        )
        step["args"][0]["visible"][i:i+num_traces] = num_traces*[True]
        steps.append(step)

    if slider == "frequency":
            value_str="Spring constant k: "
    if slider == "Delta_Q":
            value_str="&#916;Q: "
    
    sliders = [dict(
        active=0,
        currentvalue={"prefix": value_str},
        pad={"t": 50},
        steps=steps
    )]

    fig.update_yaxes(title_text="Energy",range=[-2, 30], showticklabels=True, row=1, col=1)
    fig.update_xaxes(showticklabels=False ,title_text=r"Normal coordinate Q",row=1, col=1)
    
    fig.update_xaxes(title_text="PL Intensity", row=1, col=2)
    
    fig.update_layout(
            sliders=sliders,
            showlegend=False,
            height=600, width=800)
        
    fig.show()

    