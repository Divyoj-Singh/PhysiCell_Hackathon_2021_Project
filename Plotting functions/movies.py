from pyMCDS import pyMCDS
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import seaborn as sns
sns.set_theme(context='talk', style='ticks',font='sans-serif')

folder="./../plots_NDJ/"
print("Files present in the folder")
#!ls ./output/*.xml 
#print("reading one of the files")
#mcds = pyMCDS('output00000005.xml', 'output')
L_size=250
L_size=L_size+L_size/10 
def movie_scatter_plot(t_span,centerx_over_time,centery_over_time,u_over_time,u_name="u",filename="",folder="./"):
        
    
    def plot_scatter(centerx_at_t,centery_at_t,u_at_t,t):
        # Clear the current plot figure
        plt.clf()
    
        plt.title(f"{u_name} at t = {t_span[t]:.3f}")
        plt.xlabel("x")
        plt.ylabel("y")
        
        sc=plt.scatter(centerx_at_t,centery_at_t,c=u_at_t,cmap='plasma',s=L_size/5,marker='o')
        plt.colorbar(sc)
        plt.xlim(-L_size,L_size)
        plt.ylim(-L_size,L_size)
              
        return plt
    
    def animate(t):
        
        plot_scatter(centerx_over_time[t],centery_over_time[t],u_over_time[t], t)
    
    #fig,ax=plt.figure()
    #ax.set(xlabel ='X-Axis', ylabel ='Y-Axis', xlim =(-(Lx/2)*frac, (Lx/2)*frac), ylim =(-(Ly/2)*frac, (Ly/2)*frac))
    
    anim = animation.FuncAnimation(plt.figure(),animate,interval=500,frames=np.size(u_over_time,0))
    anim.save(folder+filename+".gif",dpi=200)
    plt.close()


last_index = 50;
scatter_plot=False
times = np.zeros( last_index+1 ); 


for n in range( 0,last_index+1 ):
    
    filename='output'+"%08i"%n+'.xml'
    mcds=pyMCDS(filename,'output1')
    times[n]= mcds.get_time()
    Notch = mcds.data['discrete_cells']['notch']
    Delta = mcds.data['discrete_cells']['delta']
    Jagged = mcds.data['discrete_cells']['jagged']
    NICD = mcds.data['discrete_cells']['NICD']
        
    cx = mcds.data['discrete_cells']['position_x']
    cy = mcds.data['discrete_cells']['position_y']
    
    if scatter_plot==True:
        ## Scatter plot:
        fig, axs = plt.subplots(2,2,figsize=(16,16))        
        
        sc=axs[0,0].scatter(cx,cy,c=Notch,s=L_size/5,cmap='plasma',marker='o');
        axs[0,0].set(xlabel ='X-Axis', ylabel ='Y-Axis', xlim =(-L_size, L_size), ylim =(-L_size, L_size),title ="Notch")
        fig.colorbar(sc,ax=axs[0,0])
        
        sc=axs[0,1].scatter(cx,cy,c=Delta,s=L_size/5,cmap='plasma',marker='o');
        axs[0,1].set(xlabel ='X-Axis', ylabel ='Y-Axis', xlim = (-L_size, L_size), ylim =(-L_size, L_size),title ="Delta")
        fig.colorbar(sc,ax=axs[0,1])
        
        sc=axs[1,0].scatter(cx,cy,c=Jagged,s=L_size/5,cmap='plasma',marker='o');
        axs[1,0].set(xlabel ='X-Axis', ylabel ='Y-Axis', xlim =(-L_size, L_size), ylim =(-L_size, L_size),title ="Jagged")
        fig.colorbar(sc,ax=axs[1,0])
        
        sc=axs[1,1].scatter(cx,cy,c=NICD,s=L_size/5,cmap='plasma',marker='o');
        axs[1,1].set(xlabel ='X-Axis', ylabel ='Y-Axis', xlim =(-L_size, L_size), ylim =(-L_size, L_size),title ="NICD")
        fig.colorbar(sc,ax=axs[1,1])
        
        
        fig.suptitle("scatter plot at "+ str(mcds.get_time()))
        fig.subplots_adjust(hspace=0.5,wspace=0.5) #top=0.90, bottom=0.10, left=0.15, right=0.95, 
        fig.savefig(folder+"scatter plot at "+str(mcds.get_time())+".png")    
        plt.close()
        
    
    if n==0:
        Notch_over_time=Notch;Delta_over_time=Delta;
        Jagged_over_time=Jagged;NICD_over_time=NICD;
        centerx_over_time=cx;centery_over_time=cy;
        
    if n>0:
        Notch_over_time=np.dstack((Notch_over_time,Notch))
        Delta_over_time=np.dstack((Delta_over_time, Delta))
        Jagged_over_time=np.dstack((Jagged_over_time,Jagged ))
        NICD_over_time=np.dstack((NICD_over_time, NICD))
        centerx_over_time=np.dstack((centerx_over_time, cx))
        centery_over_time=np.dstack((centery_over_time, cy))
        #
    

Notch_over_time=np.transpose(Notch_over_time, (2, 1, 0));   
Delta_over_time=np.transpose(Delta_over_time, (2, 1, 0));   
Jagged_over_time=np.transpose(Jagged_over_time,(2, 1, 0));   
NICD_over_time=np.transpose(NICD_over_time, (2, 1, 0));  
 
centerx_over_time=np.transpose(centerx_over_time, (2, 1, 0));   
centery_over_time=np.transpose(centery_over_time, (2, 1, 0));   


movie_scatter_plot(times,centerx_over_time,centery_over_time,Notch_over_time,"Notch","notch",folder)

movie_scatter_plot(times,centerx_over_time,centery_over_time,Delta_over_time,"Delta","Delta",folder)
movie_scatter_plot(times,centerx_over_time,centery_over_time,NICD_over_time,"NICD","NICD",folder)