package org.cloudsimplus.examples.brokers;


import org.cloudbus.cloudsim.allocationpolicies.VmAllocationPolicySimple;
import org.cloudbus.cloudsim.datacenters.Datacenter;
import org.cloudbus.cloudsim.datacenters.DatacenterSimple;
import org.cloudbus.cloudsim.brokers.DatacenterBroker;
import org.cloudbus.cloudsim.cloudlets.Cloudlet;
import org.cloudbus.cloudsim.cloudlets.CloudletSimple;
import org.cloudbus.cloudsim.core.CloudSim;
//import org.cloudbus.cloudsim.datacenters.Datacenter;
import org.cloudbus.cloudsim.distributions.UniformDistr;
import org.cloudbus.cloudsim.hosts.Host;
import org.cloudbus.cloudsim.hosts.HostSimple;
//import org.cloudbus.cloudsim.power.models.PowerAware;
//import org.cloudbus.cloudsim.power.models.PowerModel;
//import org.cloudbus.cloudsim.power.models.PowerModelLinear;
import org.cloudbus.cloudsim.provisioners.PeProvisionerSimple;
import org.cloudbus.cloudsim.provisioners.ResourceProvisionerSimple;
import org.cloudbus.cloudsim.resources.Pe;
import org.cloudbus.cloudsim.resources.PeSimple;
import org.cloudbus.cloudsim.schedulers.cloudlet.CloudletSchedulerTimeShared;
import org.cloudbus.cloudsim.schedulers.vm.VmSchedulerTimeShared;
import org.cloudbus.cloudsim.utilizationmodels.UtilizationModel;
import org.cloudbus.cloudsim.utilizationmodels.UtilizationModelDynamic;
import org.cloudbus.cloudsim.utilizationmodels.UtilizationModelPlanetLab;
import org.cloudbus.cloudsim.vms.Vm;
import org.cloudbus.cloudsim.vms.VmSimple;
import org.cloudsimplus.builders.tables.CloudletsTableBuilder;
//import org.cloudsimplus.heuristics.CloudletToVmMappingHeuristic;
import org.cloudsimplus.heuristics.CloudletToVmMappingSimulatedAnnealing;
import org.cloudsimplus.heuristics.CloudletToVmMappingSolution;
//import org.cloudsimplus.heuristics.Heuristic;
//import org.cloudsimplus.util.Log;

//import ch.qos.logback.classic.Level;

import static org.cloudbus.cloudsim.brokers.DatacenterBroker.LOGGER;

import java.util.ArrayList;
//import java.util.DoubleSummaryStatistics;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;


public class DatacenterBrokerSaf {
	 private final CloudSim simulation;
	    private List<Cloudlet> cloudletList;
	    private List<Vm> vmList;
	    private CloudletToVmMappingSimulatedAnnealing heuristic;

	    /**
	     * Number of cloudlets created so far.
	     */
	    private int numberOfCreatedCloudlets = 0;
	    /**
	     * Number of VMs created so far.
	     */
	    private int numberOfCreatedVms = 0;
	    /**
	     * Number of hosts created so far.
	     */
	   // private int numberOfCreatedHosts = 0;

	    private static final int HOSTS_TO_CREATE = 1;
	    private static final int VMS_TO_CREATE = 50;
	    private static final int CLOUDLETS_TO_CREATE = 1;//number of tasks
	    public static final int   number_Of_Jobs= 50;// number of repetitive tasks
	    //comment for SA CLOUDLETS_TO_CREATE = 16 ,for DRF CLOUDLETS_TO_CREATE = 8
	    /**
		 * parameters for using planetlab
		 */
	    private static final String TRACE_FILE = "workload/planetlab/20110303/75-130-96-12_static_oxfr_ma_charter_com_irisaple_wup";
	    private static final int SCHEDULING_INTERVAL = 300;
	    /**
		 * Simulated Annealing (SA) parameters.
		 */
		public static double SA_INITIAL_TEMPERATURE = 0;
		public static final double SA_COLD_TEMPERATURE = 500;
		public static final double SA_COOLING_RATE = 0.1;
		public static final int    SA_NUMBER_OF_NEIGHBORHOOD_SEARCHES = 50;
		public static final int    number_Of_PEs= 32;
		public static final int    amount_Of_Ram = 2000 ;
		private static HashMap<Integer, Double>  fitnessOfAll = new HashMap<Integer, Double> () ;
		private static int constant = 1000;
		private static HashMap<Double, Double> distribution = new HashMap<Double,Double>(); 		
		public static final int	  number_Of_cloudlets =CLOUDLETS_TO_CREATE*number_Of_Jobs;
		public double[][] cloudletInfo = new double[number_Of_cloudlets][4];
		public static final int index_cpu_req = 0;
		public static final int index_ram_req = 1;
		public static final int index_cpu_minshare = 2;
		public static final int index_ram_minshare = 3;
		public static final double[] clusterInfo = {number_Of_PEs * HOSTS_TO_CREATE, 
												amount_Of_Ram * HOSTS_TO_CREATE};
		public int[] dominant = new int[number_Of_cloudlets];
		public double[][] shares_Cloudlet = new double[number_Of_cloudlets][2];
		public double[][] usage_Cloudlet = new double[number_Of_cloudlets][2];
		public ArrayList<Boolean> isNeedy= new ArrayList<Boolean>();
		public org.cloudbus.cloudsim.brokers.DatacenterBrokerHeuristic broker0;
		public static Map<Cloudlet, Vm> cloudletVmMap1 = new HashMap<Cloudlet, Vm>();
		double[][] frees = new double[VMS_TO_CREATE][2];
		double cost;
		double fitness;	
		List<Host> hostList = new ArrayList<>();
		int pasvand=-1;
		//**************************************************************************
		
		public static void main(String[] args) {
	        new DatacenterBrokerSaf();
	    }
		
		private static  double CalculateTemperature() {
			double sum = 0;
	    	  
	    	  HashMap< Double, Double> distribution = calculatePdf(fitnessOfAll);
	    	  Iterator<Map.Entry< Double, Double>> it = distribution.entrySet().iterator();
	    	    while (it.hasNext()) {
	    	    	HashMap.Entry< Double, Double> pair = (HashMap.Entry< Double, Double>)it.next();
	    	    	double logOfPi= Math.log10((double) pair.getValue());
	    	    	double Temperature =  -1 *constant  * logOfPi * ((double) pair.getValue());
	    	    	sum += Temperature ;
	    	    	
	    	    }
	    	  return sum;
			
		}
		
		//*********************************************************************
		static HashMap<Double, Double> calculatePdf(HashMap<Integer, Double> fitnessofall){
	    	  double totalNumOfFitnesses = 0;
	    	  Iterator<HashMap.Entry<Integer,Double>> itr = fitnessofall.entrySet().iterator();
	    	    while (itr.hasNext()) {
	    	    	HashMap.Entry<Integer,Double> pair = (HashMap.Entry<Integer,Double>)itr.next();
	    	    	double key = (double) pair.getValue();
	    		  if (distribution.containsKey(key))
	    		  {
	                double count = distribution.get(key);
	                count++;
	                distribution.put( key, count );
	    		  }else
	            {
	            	distribution.put( key, (double) 1 );
	            }
	    		  
	    	  }
	    	    
	    	  Iterator<HashMap.Entry<Double,Double>> it = distribution.entrySet().iterator();
	    	    while (it.hasNext()) {
	    	    	HashMap.Entry<Double,Double> pair = (HashMap.Entry<Double,Double>)it.next();
	    	    	totalNumOfFitnesses += (double)pair.getValue();
	    	    }
	    	    it = distribution.entrySet().iterator();
	    	    HashMap<Double, Double> finalDistribution =new HashMap<Double, Double>();
	    	    while (it.hasNext()) {
	    	    	HashMap.Entry<Double,Double> pair = (HashMap.Entry<Double,Double>)it.next();
	    	    	finalDistribution.put((double) pair.getKey(), (double)pair.getValue() / totalNumOfFitnesses );
	    	    }
	    	    
		return finalDistribution;
	      }
		//****************************************************************************
		
		 void calculateClusterAndFairRatios(double[][] cloudletInfo,
			        float weight, double[] clusterInfo) {
			      
			      int index_ram = index_ram_req;
			      int index_cpu = index_cpu_req;
                  for(int i=0; i< CLOUDLETS_TO_CREATE; i++) {
			      shares_Cloudlet[i][index_ram] =
			          ((double) cloudletInfo[i][index_ram_req]) /
			          clusterInfo[index_ram];
			      shares_Cloudlet[i][index_cpu] =
			          ((double) cloudletInfo[i][index_cpu_req]) /
			          clusterInfo[index_cpu];
			      dominant [i] =
			    		  shares_Cloudlet[i][index_cpu] > shares_Cloudlet[i][index_ram] ?
			        		  index_cpu : index_ram;

			    		  shares_Cloudlet[i][index_ram] /= weight;
			    		  shares_Cloudlet[i][index_cpu] /= weight;

			     }
			    }  

		
		
		//****************************************************************************
		 
		 void UpdateNeedy (ArrayList<Cloudlet> cloudlets) {
			 for (Cloudlet C: cloudlets) {
				 if( usage_Cloudlet[(int) C.getId()][dominant[(int) C.getId()]] <
						 cloudletInfo[(int) C.getId()][dominant[(int) C.getId()]+2]) {
					 
					 isNeedy.add((int) C.getId(), true);
				 }else 
				 {isNeedy.add((int) C.getId(), false);
				 
				 }
			}
			 
			  
		 }
		  //**************************************************************************
		   void DRfScheduling(ArrayList<Cloudlet> cloudlets) {	
			  ArrayList<Cloudlet> temp = new ArrayList<Cloudlet>();
			  Map<Integer, Double> unSortedMap = new HashMap<Integer, Double>();
			  
			  for (Cloudlet c:cloudlets) {  
				  if (isNeedy.get((int) c.getId()).equals(true) ){
					  temp.add((int) c.getId(),c);  // fill temp with needy cloudlets
				  }
				  
			  }
				
			  for (Cloudlet c:temp) {
				  //fill sorted but not sorted yet  (cloudlet id, drf shares)
				  unSortedMap.put((int)c.getId(), calculateDRFShares(c,1));
			  }
			  
			    //start ascending sort 
			   LinkedHashMap<Integer, Double> sortedMap = new LinkedHashMap<>();
		       unSortedMap.entrySet().stream().sorted(Map.Entry.comparingByValue())
		                .forEachOrdered(x -> sortedMap.put(x.getKey(), x.getValue()));
		       
			  //bind to VM pick with lowest dominant share
		          for (Map.Entry<Integer, Double> entry : sortedMap.entrySet()) {
		    	 //bind the cloudlets to the vms. This way, the broker
		           // will submit the bound cloudlets only to the specific VM
		          	  
		             Cloudlet cloudlet = cloudletList.get(entry.getKey());
		             Vm bestVm= null;
		             if (!broker0.getCloudletFinishedList().contains(cloudlet)) {
		               //filter vm that fits cloudlets requirement.
		            	 bestVm=bestFitCloudletToVm(cloudlet)	;   
		            	 bindCloudletToVm(cloudlet, bestVm);
		            	 AddTOArray(cloudlet);		            	 
				        }
		                       
		       	       }
		         			            
			}
		   //**********************************************************************
		   
		   //**********************************************************************
		  private void bindCloudletToVm(Cloudlet cloudlet, Vm bestVm) {
			  cloudletVmMap1.put(cloudlet, bestVm);
			  updatefreespace(cloudlet, bestVm);
		}
		  //**************************************************************************
		private void updatefreespace(Cloudlet cloudlet, Vm bestVm) {
			
	    	 frees[(int)bestVm.getId()][0]-= (cloudlet.getNumberOfPes());
	    	 frees[(int)bestVm.getId()][1]-= (cloudlet.getAmountOfRam());   	
	    	    		 
	    	    	}
			
		//***********************************************************************
		   private Vm bestFitCloudletToVm(final Cloudlet cloudlet) {
			   ArrayList<Vm> tempVm = new ArrayList<Vm>();
			  
				for(Vm vm:vmList) 
				{
						if( frees[(int)vm.getId()][0] >= cloudlet.getNumberOfPes() && 
							frees[(int)vm.getId()][1] >= cloudlet.getAmountOfRam()) 
						{
						tempVm.add(vm);
						}
						
				}
				
				double minPesOfVms= number_Of_PEs;
				Vm minVm=null;
				//to resolve nullpointerexception
				/*if (tempVm.size()!=0) {				
				 minVm= tempVm.get( getRandomNumberOfPes(tempVm.size()-1));	}
				else {
				 minVm= vmList.get(0);	
				}*/
					for(Vm vm:tempVm)
				{  double VmnumberPes= vm.getNumberOfPes();
					 if( VmnumberPes <= minPesOfVms) 
					 {
						 minPesOfVms=vm.getNumberOfPes(); 
						 minVm= vm;
					 }
						 
				}	
				return minVm;
		     }
		   
		 //***********************************************************************
		   @SuppressWarnings("unused")
		   private Vm firstFitCloudletToVm(final Cloudlet cloudlet) {
			   ArrayList<Vm> tempVm = new ArrayList<Vm>();	
			   pasvand = getRandomAmountOfRam(VMS_TO_CREATE);
			  while (frees[pasvand][0] < cloudlet.getNumberOfPes() || 
					     frees[pasvand][1] < cloudlet.getAmountOfRam()) {
				  pasvand = getRandomAmountOfRam(VMS_TO_CREATE);
			  }			  
				for(Vm vm:vmList) 
				{
						if( (int)vm.getId()==pasvand) 
						{ tempVm.add(vm);}			 			
										
				}			
				
				Vm minVm=null;
					for(Vm vm:tempVm)
				         {   minVm= vm;  }	 				 
				return minVm;
		     }
		   
		
		//****************************************************************************
		  double calculateDRFShares(Cloudlet c,float weight) {
			        // drfFairshare = request(dominant)/minshare(dominant)*weight
			  double drfFairshare =  0;
			  drfFairshare  = 
					  cloudletInfo[(int)c.getId()][dominant[(int)c.getId()]] / 
					  cloudletInfo[(int)c.getId()][dominant[(int)c.getId()]+2] * weight;  
			  return drfFairshare;
		  }
			      
		//****************************************************************************
		double calculateOurFairness(Cloudlet c, float weight) {
		       //ourfairness=usage(dominant)/clustercapacity(dominant)*weight
		      
		    double ourFairness =  0;
			ourFairness =
					 usage_Cloudlet[(int)c.getId()][dominant[(int)c.getId()]] /
					clusterInfo[dominant[(int)c.getId()]] * weight;   

		   return ourFairness;
	   }
		
	    //**************************************************************************

		 float calculateFitnessOfCloudlet(Cloudlet C) {
			 double[] clusterCapacity=getClusterCapacity() ;
			 	  double CoreCapacity =clusterCapacity[0];
			 	  double RamCapacity =clusterCapacity[1];
			 	  //fitness= request * available capacity
			 	  float fitness=(float) (( C.getAmountOfRam()* RamCapacity) +
			     		( C.getNumberOfPes()* CoreCapacity)) ;
			      return fitness;
			    }
		    //**************************************************************************


		 double[] getClusterCapacity() {
			 double[] clusterCapacity=new double[2] ;
			 	  for(Vm vm :vmList) {
			 		//clusterCapacity = summation of first column of frees 
			 	clusterCapacity[index_cpu_req]+= frees[(int)vm.getId()][0]; //cpu
			 	clusterCapacity[index_ram_req]+= frees[(int)vm.getId()][1]; //ram
			 	}
				return clusterCapacity;
			    }
		    //**************************************************************************

		 void AddTOArray(Cloudlet c) {
	    	  //if c does'nt exist
			 double fitness =calculateFitnessOfCloudlet(c);
			 
			 if(!fitnessOfAll.containsKey(c.getId())) //cloudlet doesnt exist
			 {
	 			fitnessOfAll.put((int)c.getId(),fitness);
			
			 }else if(!broker0.getCloudletFinishedList().contains(c)) // cloudlet exists so update 
			 {
				 fitnessOfAll.replace((int)c.getId(),fitness);
				 
			 }else if(broker0.getCloudletFinishedList().contains(c)) 
			 {
				 fitnessOfAll.remove(c.getId(),fitness);   //remove cloudlet if its finished
			 }
	 	
	 	}
		 
		    //**************************************************************************
		 
		public DatacenterBrokerSaf() {
			//initialize Free spaces of vcores and ram
	        for (int i=0; i < VMS_TO_CREATE ;++i) {
				 frees [i][0]=number_Of_PEs;
				 frees [i][1]=amount_Of_Ram;
			 }	       
	        System.out.println("Starting " + getClass().getSimpleName());
	        this.vmList = new ArrayList<>();
	        this.cloudletList = new ArrayList<>();

	        simulation = new CloudSim();
	        @SuppressWarnings("unused")
	        Datacenter datacenter0 = createDatacenter();
	        broker0 = createBroker();      
	        createAndSubmitVms(broker0);
	        createAndSubmitCloudlets(broker0);	
	        final long startTime = System.currentTimeMillis();
	        calculateClusterAndFairRatios(cloudletInfo, 1, clusterInfo);
	        Scheduling();	        
	     // simulation.addOnClockTickListener(this::dynamicCloudletArrival);
	        broker0.setVmMapper(this::Function);
	        simulation.start();	 	        
	        List<Cloudlet> finishedCloudlets = broker0.getCloudletFinishedList();
	        new CloudletsTableBuilder(finishedCloudlets).build();
	        printCloudlets();	        
	        final double finishtime=( (System.currentTimeMillis() - startTime)/1000.0);
	        LOGGER.info(
	                " finished the solution find for mapping Cloudlets to Vm's in\n {} seconds with a solution cost of {}",
	                 finishtime);	        
	     	print(broker0);
	     	
	    }
	//******************************************************************************   
	/*	private void dynamicCloudletArrival(EventInfo e) {
		        if((int)e.getTime() != 35) {
		            return;
		        }
				
			}*/
		//**************************************************************************
	private Vm Function (Cloudlet cloudlet1) {
	 return cloudletVmMap1.get(cloudlet1);
        		}
       //****************************************************************************
		
		private org.cloudbus.cloudsim.brokers.DatacenterBrokerHeuristic createBroker() {
			createSimulatedAnnealingHeuristic();
			org.cloudbus.cloudsim.brokers.DatacenterBrokerHeuristic broker0 = new org.cloudbus.cloudsim.brokers.DatacenterBrokerHeuristic(simulation);
			broker0.setHeuristic(heuristic);
			return broker0;
		}
       
		private void createAndSubmitCloudlets(org.cloudbus.cloudsim.brokers.DatacenterBrokerHeuristic broker0) {
			int k=0;
			double clocktime = 0;				
			for(int i = 0; i < CLOUDLETS_TO_CREATE; i++){
			//	int cpu_minshare =getRandomNumberOfPes(number_Of_PEs / 2); 
			//	int ram_minshare =getRandomAmountOfRam(amount_Of_Ram / 2); 
				int cpu_temp1 =number_Of_PEs/2;
				int cpu_temp2 =2;
						//getRandomNumberOfPes(number_Of_PEs);
				int ram_temp1 =amount_Of_Ram;
				int ram_temp2 =100;
						//getRandomAmountOfRam(amount_Of_Ram-513)+512;
				clocktime += 10;
						//getRandomAmountOfRam(180);
				//cpu_temp++;
				//ram_temp=ram_temp+10;				
			
				for(int j=0; j<number_Of_Jobs;j++) {	
				if (j%2==0) {		
					if(cpu_temp2+j>(number_Of_PEs/2)) {cpu_temp2-=((number_Of_PEs/2)-1);}
					if(ram_temp2+(j*51)>amount_Of_Ram) {ram_temp2-=(amount_Of_Ram-100);}
			   cloudletList.add(createCloudlet(broker0, cpu_temp2+j,ram_temp2+(j*51),clocktime)) ;
			   cloudletInfo[k][index_cpu_req]= cpu_temp2+j;
			   cloudletInfo[k][index_ram_req]= ram_temp2+(j*51);
			   cloudletInfo[k][index_cpu_minshare]= (cpu_temp2+j)/2;				   
			   cloudletInfo[k][index_ram_minshare]= (ram_temp2+(j*51))/2;					   
			   usage_Cloudlet[k][0]= cpu_temp2+j;
   			   usage_Cloudlet[k][1]= ram_temp2+(j*51);
			   k++;			  
			   }else {
				   if (cpu_temp1-j<1) {cpu_temp1+=((number_Of_PEs/2)-1);}
				   if (ram_temp1-(j*51)<10) {ram_temp1+=(amount_Of_Ram-100);}
				   cloudletList.add(createCloudlet(broker0, cpu_temp1-j,ram_temp1-(j*51),clocktime)) ;
				   cloudletInfo[k][index_cpu_req]= cpu_temp1-j;
				   cloudletInfo[k][index_ram_req]= ram_temp1-(j*51);
				   cloudletInfo[k][index_cpu_minshare]= (cpu_temp1-j)/2;				   
				   cloudletInfo[k][index_ram_minshare]= (ram_temp1-(j*51))/2;					   
				   usage_Cloudlet[k][0]= cpu_temp1-j;
	   			   usage_Cloudlet[k][1]= ram_temp1-(j*51);
				   k++; 				   
			   }
				
				}		
		}	
			broker0.submitCloudletList(cloudletList);
		}
	    //**************************************************************************
	
		
		private void createAndSubmitVms(org.cloudbus.cloudsim.brokers.DatacenterBrokerHeuristic broker0) {
			vmList = new ArrayList<>(VMS_TO_CREATE);
			for(int i = 0; i < VMS_TO_CREATE; i++){
			    vmList.add(createVm(broker0, number_Of_PEs));   //number of pes
			}
			broker0.submitVmList(vmList);
		}
	    //**************************************************************************

		private void createSimulatedAnnealingHeuristic() {
			SA_INITIAL_TEMPERATURE=903.0899869919433;
			heuristic =
			        new CloudletToVmMappingSimulatedAnnealing(SA_INITIAL_TEMPERATURE, new UniformDistr(0, 1));
			heuristic.setColdTemperature(SA_COLD_TEMPERATURE);
			heuristic.setCoolingRate(SA_COOLING_RATE);
			heuristic.setNeighborhoodSearchesByIteration(SA_NUMBER_OF_NEIGHBORHOOD_SEARCHES);
		}
	    //**************************************************************************
		private void printCloudlets() {
			System.out.println("**********Cluster Capacity:***********");
			System.out.printf("ClusterCapacity: <%fCpu,%fRam>",clusterInfo[0],clusterInfo[1]);
			System.out.println("\n\n**********Cloudlet Requests:***********");
			for(Cloudlet C: cloudletList) 
			{
				System.out.printf("Cloudlet Id: %f |Request: <%fCpu,%f Ram>| minshare:"
						+ "<%fCpu,%f Ram>\n",(float) C.getId(),
						(float)cloudletInfo[(int) C.getId()][0],
						(float)cloudletInfo[(int) C.getId()][1],
						(float)cloudletInfo[(int) C.getId()][2],
						(float)cloudletInfo[(int) C.getId()][3]);
								
			}
			System.out.println("***********************************************************************");
		}
		//*****************************************************************************
		private void print(org.cloudbus.cloudsim.brokers.DatacenterBrokerHeuristic broker0) {
							
			 CloudletToVmMappingSolution FairnessSolution = new CloudletToVmMappingSolution(heuristic, cloudletVmMap1);
		
			 System.out.printf(" (cost: %.2f fitness: %.6f)\n",
					 FairnessSolution.getCost(),  FairnessSolution.getFitness());
								
			System.out.println("Simulated Annealing Parameters");
			System.out.printf("\tInitial Temperature: %.2f", SA_INITIAL_TEMPERATURE);
			System.out.printf(" Cooling Rate: %.4f", SA_COOLING_RATE);
			System.out.printf(" Cold Temperature: %.6f", SA_COLD_TEMPERATURE);
			//System.out.printf(" Number of neighborhood searches by iteration: %d\n", SA_NUMBER_OF_NEIGHBORHOOD_SEARCHES);
	        System.out.println(getClass().getSimpleName() + " finished!");
		}
	    //**************************************************************************


		/**
	     * Randomly gets a number of PEs (CPU cores).
	     *
	     * @param maxPesNumber the maximum value to get a random number of PEs
	     * @return the randomly generated PEs number
	     */
		@SuppressWarnings("unused")
	    private int getRandomNumberOfPes(int maxPesNumber) {
	        return heuristic.getRandomValue(maxPesNumber)+1;
	    }
	    //**************************************************************************
	    /**
	     * Randomly gets a amount of ram (CPU cores).
	     *
	     * @param maxRamAmount the maximum value to get a random amount of ram
	     * @return the randomly generated Ram amount
	     */		
	    private int getRandomAmountOfRam(int maxRamAmount) {
	        return heuristic.getRandomValue(maxRamAmount);
	    }
	    //**************************************************************************
  	// to create datacenter
	        private DatacenterSimple createDatacenter() {
	        	
	        for(int i = 0; i < HOSTS_TO_CREATE; i++) {
	            hostList.add(createHost());
	        }

	        return new DatacenterSimple(simulation, hostList, new VmAllocationPolicySimple());
	    }
	    //**************************************************************************


	    public Host createHost() {
	        long mips = 1000*VMS_TO_CREATE; // capacity of each CPU core (in Million Instructions per Second)
	       // int  ram = 8192; // host memory (Megabyte)
	        long storage = 1000000* VMS_TO_CREATE; // host storage
	        long bw = 10000*VMS_TO_CREATE;	        
	        List<Pe> peList = new ArrayList<>();
	        /*Creates the Host's CPU cores and defines the provisioner
	        used to allocate each core for requesting VMs.*/
	        for(int i = 0; i < 8; i++)
	            peList.add(new PeSimple(mips, new PeProvisionerSimple()));

	       return new HostSimple(amount_Of_Ram*VMS_TO_CREATE, bw, storage, peList)	    		  
	           .setRamProvisioner(new ResourceProvisionerSimple())
	           .setBwProvisioner(new ResourceProvisionerSimple())
	            .setVmScheduler(new VmSchedulerTimeShared());
	    }
	    //**************************************************************************

	    private Vm createVm(DatacenterBroker broker, int pesNumber) {
	        long mips = 100;
	        long   storage = 1000000; // vm image size (Megabyte)
	       // int    ram = 512; // vm memory (Megabyte)
	        long   bw = 10000; // vm bandwidth

	        return new VmSimple(numberOfCreatedVms++, mips, pesNumber)
	            .setRam(amount_Of_Ram).setBw(bw).setSize(storage)
	            .setCloudletScheduler(new CloudletSchedulerTimeShared());

	    }
	    //**************************************************************************


	    private Cloudlet createCloudlet(DatacenterBroker broker, int numberOfPes , 
	    								int amountOfRam, double clocktime) {	    	
	    	int CLOUDLET_LENGTH = 2000;
	    	//to use planetlab entries
	    	final UtilizationModel utilizationCpu = UtilizationModelPlanetLab.getInstance(TRACE_FILE, SCHEDULING_INTERVAL);
	           Cloudlet cloudlet =
	             new CloudletSimple(numberOfCreatedCloudlets++, CLOUDLET_LENGTH,numberOfPes,amountOfRam)
	                    .setFileSize(4096)
	                    .setOutputSize(4096)
	                    .setUtilizationModelCpu(utilizationCpu)
	                    .setUtilizationModelBw(new UtilizationModelDynamic(0.99))
	                    .setUtilizationModelRam(new UtilizationModelDynamic(0.99))
	                    ;
	           cloudlet.setExecStartTime(clocktime);
	            
	    	return cloudlet;
	    	/* long length = 400000; //in Million Structions (MI)
	        long fileSize = 300; //Size (in bytes) before execution
	        long outputSize = 300; //Size (in bytes) after execution

	        //Defines how CPU, RAM and Bandwidth resources are used
	        //Sets the same utilization model for all these resources.
	        UtilizationModel utilization = new UtilizationModelFull();

	        return new CloudletSimple(numberOfCreatedCloudlets++, length, numberOfPes, amountOfRam)
	        .setFileSize(fileSize)
	        .setOutputSize(outputSize)
	        .setUtilizationModel(utilization);*/
	    }
	    //**************************************************************************

	    private double Scheduling() {
	    	     	  	
	           //determining cloudlets that are needy or not  
	           UpdateNeedy((ArrayList<Cloudlet>) cloudletList);
	           
	            if (isNeedy.contains(true))
	            {
	            	//drf phase1(select cloudlet with lowest minshare
	            	DRfScheduling((ArrayList<Cloudlet>) cloudletList);
	                    	
	            }
	            else {	
	            	//simulated annealing & drf phase2(select cloudlet with lowest drf shares)
	            	SimulatedAnnealingScheduling((ArrayList<Cloudlet>) cloudletList);
	            }
	            
	           
			return 0;
	    }
	  
		//**************************************************************************

	    private void SimulatedAnnealingScheduling(ArrayList<Cloudlet> cloudletList) {
	        	
	    	 Map<Integer, Double> unSortedMap = new HashMap<Integer, Double>() ;
	    	for (Cloudlet c:cloudletList) {
				  //fill sorted but not sorted yet  (cloudlet id, drf shares or ourfairness) 
				  unSortedMap .put((int)c.getId(), calculateOurFairness(c,1));
			  }
			  
			/*    //start ascending sort 
			   LinkedHashMap<Integer, Double> sortedMap = new LinkedHashMap<>();
		       unSortedMap.entrySet().stream().sorted(Map.Entry.comparingByValue())
		                .forEachOrdered(x -> sortedMap.put(x.getKey(), x.getValue())); */
		       
		    // uncomment for test calculate temperature => double s=0;
		        // start Simulated Annealing
		       //uncomment for SA
             createSimulatedAnnealingHeuristic(); 
		       double T = SA_INITIAL_TEMPERATURE;  
		       Iterator<Map.Entry<Integer, Double>> it = unSortedMap.entrySet().iterator();
		       while (it.hasNext()) {
		    	   Map.Entry<Integer, Double> pair = (Map.Entry<Integer, Double>)it.next();		        		
		    	   int  key1= (int) pair.getKey();
		    	   double  value1 = (double) pair.getValue();	   	
		    	   pair = (Map.Entry<Integer, Double>)it.next();
		    	   int  key2= (int) pair.getKey();
		    	   double  value2 = (double) pair.getValue();
		    	   	 Cloudlet cloudlet1 = cloudletList.get(key1);
		        	 Cloudlet cloudlet2 = cloudletList.get(key2);
		        	 Cloudlet cloudlet = cloudlet1;
		        	 Cloudlet cloudleth = cloudlet2;
		        	 Vm bestVm= null;
		        	 Vm bestVmh= null;
		        	 double p= Math.exp( ( (value2) - (value1) ) / T  );
		    	      	  if(p>Math.random()) {
		    	      		cloudlet = cloudlet2;
		    	      		cloudleth = cloudlet1;
		    	      			      	 }
		    	   	if (!broker0.getCloudletFinishedList().contains(cloudlet)) {
		    	   		//pasvand = getRandomAmountOfRam(VMS_TO_CREATE);
		    	      	bestVm=bestFitCloudletToVm(cloudlet);
		    	      	bindCloudletToVm(cloudlet, bestVm);
		    	      	bestVmh=bestFitCloudletToVm(cloudleth);
		    	      	bindCloudletToVm(cloudleth, bestVmh);
		    	      	AddTOArray(cloudlet);
		    	      	AddTOArray(cloudleth);
		        	    if (T > SA_COLD_TEMPERATURE) {		        	    	
		        	    T *= (1- SA_COOLING_RATE) ;
		        	    }else {T=CalculateTemperature();}
		        	    	System.out.printf("\n*********** Temperature: %.2f", T);
		       	     
		         }}
	    	// uncomment for test calculate temperature => s=CalculateTemperature(); 		          
		      // uncomment for DRF
		     
		       	 //bind the cloudlets to the vms. This way, the broker
		           // will submit the bound cloudlets only to the specific VM	    	
	 	/*    Iterator<Map.Entry<Integer, Double>> it = unSortedMap.entrySet().iterator();
			       while (it.hasNext()) {
			    	   Map.Entry<Integer, Double> pair = (Map.Entry<Integer, Double>)it.next();		        		
			    	   int  key1= (int) pair.getKey();
			    	   double  value1 = (double) pair.getValue();	   	
			    	   pair = (Map.Entry<Integer, Double>)it.next();
			    	   int  key2= (int) pair.getKey();
			    	   double  value2 = (double) pair.getValue();
			    	   	 Cloudlet cloudlet1 = cloudletList.get(key1);
			        	 Cloudlet cloudlet2 = cloudletList.get(key2);
			        	 Cloudlet cloudlet = cloudlet1;
			        	 Cloudlet cloudleth = cloudlet2;
			        	 Vm bestVm= null;
			        	 Vm bestVmh= null;			        	
			    	        if((value2) < (value1)) {
			    	      		cloudlet = cloudlet2;
			    	      		cloudleth = cloudlet1;
			    	      			      	 }
			    	      	if (!broker0.getCloudletFinishedList().contains(cloudlet)) {			    	      		
				    	      	bestVm=bestFitCloudletToVm(cloudlet);
				    	      	bindCloudletToVm(cloudlet, bestVm);				    	      	
				    	      	bestVmh=bestFitCloudletToVm(cloudleth);
				    	      	bindCloudletToVm(cloudleth, bestVmh);
				    	      	AddTOArray(cloudlet);
				    	      	AddTOArray(cloudleth);
			    	      	}
		             		                       
		       	       }    */
	    }		
		
}
