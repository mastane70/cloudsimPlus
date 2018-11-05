package org.cloudsimplus.examples.brokers;


import org.cloudbus.cloudsim.allocationpolicies.VmAllocationPolicySimple;
import org.cloudbus.cloudsim.brokers.DatacenterBroker;
import org.cloudbus.cloudsim.brokers.DatacenterBrokerHeuristic;
import org.cloudbus.cloudsim.cloudlets.Cloudlet;
import org.cloudbus.cloudsim.cloudlets.CloudletSimple;
import org.cloudbus.cloudsim.core.CloudSim;
import org.cloudbus.cloudsim.datacenters.Datacenter;
import org.cloudbus.cloudsim.datacenters.DatacenterSimple;
import org.cloudbus.cloudsim.distributions.ContinuousDistribution;
import org.cloudbus.cloudsim.distributions.UniformDistr;
import org.cloudbus.cloudsim.hosts.Host;
import org.cloudbus.cloudsim.hosts.HostSimple;
import org.cloudbus.cloudsim.provisioners.PeProvisionerSimple;
import org.cloudbus.cloudsim.provisioners.ResourceProvisionerSimple;
import org.cloudbus.cloudsim.resources.Pe;
import org.cloudbus.cloudsim.resources.PeSimple;
import org.cloudbus.cloudsim.schedulers.cloudlet.CloudletSchedulerTimeShared;
import org.cloudbus.cloudsim.schedulers.vm.VmSchedulerTimeShared;
import org.cloudbus.cloudsim.utilizationmodels.UtilizationModel;
import org.cloudbus.cloudsim.utilizationmodels.UtilizationModelFull;
import org.cloudbus.cloudsim.vms.Vm;
import org.cloudbus.cloudsim.vms.VmSimple;
import org.cloudsimplus.builders.tables.CloudletsTableBuilder;
import org.cloudsimplus.heuristics.CloudletToVmMappingSimulatedAnnealing;
import org.cloudsimplus.heuristics.CloudletToVmMappingSolution;

import ch.qos.logback.classic.pattern.Util;

import java.nio.file.attribute.UserDefinedFileAttributeView;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Optional;
import java.util.SortedSet;
import java.util.TreeSet;

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
	    private int numberOfCreatedHosts = 0;

	    private static final int HOSTS_TO_CREATE = 1;
	    private static final int VMS_TO_CREATE = 8;
	    private static final int CLOUDLETS_TO_CREATE = 2;
	    
		/**
		 * Simulated Annealing (SA) parameters.
		 */
		public static double SA_INITIAL_TEMPERATURE = 0;
		public static final double SA_COLD_TEMPERATURE = 0.0001;
		public static final double SA_COOLING_RATE = 0.003;
		public static final int    SA_NUMBER_OF_NEIGHBORHOOD_SEARCHES = 50;
		public static final int    number_Of_PEs= 8;
		public static final int    amount_Of_Ram = 1000 ;
		private static HashMap<Integer, Double>  fitnessOfAll = new HashMap<Integer, Double> () ;
		private static int constant = 1000;
		//private static HashMap<Float,Float> usedVmResources = new HashMap<Float,Float>();
		private static HashMap<Integer, Double> distribution = new HashMap<Integer,Double>(); 
		public static final int    number_Of_Jobs= 4;
		public static final int		number_Of_cloudlets =CLOUDLETS_TO_CREATE*number_Of_Jobs;
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
	    //**************************************************************************
		
		public static void main(String[] args) {
	        new DatacenterBrokerSaf();
	    }
		private static  double CalculateTemperature() {
			// TODO Auto-generated method stub
			 double sum = 0;
	    	  
	    	  HashMap<Integer, Double> distribution = calculatePdf(fitnessOfAll);
	    	  Iterator it = distribution.entrySet().iterator();
	    	    while (it.hasNext()) {
	    	    	HashMap.Entry pair = (HashMap.Entry)it.next();
	    	    	double logOfPi= Math.log10((double) pair.getValue());
	    	    	double Temperature =  -1 *constant  * logOfPi * ((double) pair.getValue());
	    	    	sum += Temperature ;
	    	    	
	    	    }
	    	  return sum;
			
		}
		
		//*********************************************************************
		static HashMap<Integer, Double> calculatePdf(HashMap<Integer, Double> fitnessofall){
	    	  double totalNumOfFitnesses = 0;
	    	  Iterator itr = distribution.entrySet().iterator();
	    	    while (itr.hasNext()) {
	    	    	HashMap.Entry pair = (HashMap.Entry)itr.next();
	    	    	double key = (double) pair.getValue();
	    		  if (distribution.containsKey(key))
	    		  {
	                double count = distribution.get(key);
	                count++;
	                distribution.put((int) key, count );
	    		  }else
	            {
	            	distribution.put((int) key, (double) 1 );
	            }
	    		  
	    	  }
	    	    
	    	  Iterator it = distribution.entrySet().iterator();
	    	    while (it.hasNext()) {
	    	    	HashMap.Entry pair = (HashMap.Entry)it.next();
	    	    	totalNumOfFitnesses += (double)pair.getValue();
	    	    }
	    	    it = distribution.entrySet().iterator();
	    	    HashMap<Integer, Double> finalDistribution =new HashMap<Integer, Double>();
	    	    while (it.hasNext()) {
	    	    	HashMap.Entry pair = (HashMap.Entry)it.next();
	    	    	finalDistribution.put((Integer) pair.getKey(), (double)pair.getValue() / totalNumOfFitnesses );
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
				// TODO Auto-generated method stub
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
		       
		       int i = 0;
		       for (Map.Entry<Integer, Double> entry : sortedMap.entrySet()) {
		    	 //bind the cloudlets to the vms. This way, the broker
		           // will submit the bound cloudlets only to the specific VM
		          
		    	   
		             Cloudlet cloudlet = cloudletList.get(entry.getKey());
		             
		             if (!broker0.getCloudletFinishedList().contains(cloudlet)) {
		            	 			//filter vm that fits cloudlets requirement.
		            			ArrayList<Vm> tempVm = new ArrayList<Vm>();
		            			for(Vm vm:vmList) 
		            			{
		            					double freeVmCpu= vm.getNumberOfPes()*(1-vm.getCpuPercentUsage()*0.01);
		            					double freeVmRam=vm.getRam().getAvailableResource();
		            					if(freeVmCpu >= cloudlet.getNumberOfPes() && 
		            						freeVmRam >= cloudlet.getAmountOfRam()) 
		            					{
		            					tempVm.add(vm);
		            					}
		            			}
		            			double minPesOfVms=tempVm.get(0).getNumberOfPes();
		            			Vm minVm= null;
		            			for(Vm vm:tempVm)
		            			{  double VmnumberPes= vm.getNumberOfPes();
		            				 if( VmnumberPes <= minPesOfVms) 
		            				 {
		            					 minPesOfVms=vm.getNumberOfPes(); 
		            					 minVm= vm;
		            				 }
		            					 
		            			}	
		           cloudlet.setVm(minVm);	
		           AddTOArray(cloudlet);		
		            			
	/*cloudlet.getBroker().getVmCreatedList().stream().filter(vm -> vm.getNumberOfPes() >= cloudlet.getNumberOfPes())
		     	                .min(Comparator.comparingLong(Vm::getNumberOfPes)).orElse(vmList.get(7));*/
		            	   
				        }
		    	 
		       
		       }		 
			            
			}
		  //************************************************************************
		  
		//****************************************************************************
		  double calculateDRFShares(Cloudlet c,float weight) {
			        // drfFairshare = usage(dominant)/minshare(dominant)*weight
			  double drfFairshare =  0;
			  drfFairshare  = 
					  usage_Cloudlet[(int)c.getId()][dominant[(int)c.getId()]] / 
					  cloudletInfo[(int)c.getId()][dominant[(int)c.getId()]+2] * weight;  
			  return drfFairshare;
		  }
			      
		//****************************************************************************
		double calculateOurFairness(Cloudlet c, float weight) {
		       
		      
		    double ourFairness =  0;
			ourFairness =
					 usage_Cloudlet[(int)c.getId()][dominant[(int)c.getId()]] /
					clusterInfo[dominant[(int)c.getId()]] * weight;   

		   return ourFairness;
	   }
		
	    //**************************************************************************

		 float calculateFitnessOfCloudlet(Cloudlet C) {
			 double[] clusterCapacity=getClusterCalpacity() ;
			 	  double vmCoreCapacity =clusterCapacity[0];
			 	  double vmRamCapacity =clusterCapacity[1];
			 	   //edit
			      float fitness=(float) (( C.getAmountOfRam()* vmRamCapacity) +
			     		( C.getNumberOfPes()* vmCoreCapacity)) ;
			      return fitness;
			    }
		    //**************************************************************************


		 double[] getClusterCalpacity() {
			 double[] clusterUsage=new double[2] ;
			 double[] clusterCapacity=new double[2] ;
			 	  for(Vm vm :broker0.getVmCreatedList()) {
			 		//usage percent of vm * total pes of vm= number of used pe's
			 		 clusterUsage[0] += vm.getCpuPercentUsage()*vm.getNumberOfPes(); //CPU
			 				clusterUsage[1] += vm.getRam().getAvailableResource(); //Ram
			 		
			 	  }
			 	 clusterCapacity[index_cpu_req]= clusterInfo[index_cpu_req] - clusterUsage[index_cpu_req];
			 	clusterCapacity[index_ram_req]= clusterInfo[index_ram_req] - clusterUsage[index_ram_req];
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
				 fitnessOfAll.remove(c.getId());   //remove cloudlet if its finished
			 }
	 	
	 	}
		 
		    //**************************************************************************
		 
		public DatacenterBrokerSaf() {
	        /*Enables just some level of log messages.
	          Make sure to import org.cloudsimplus.util.Log;*/
	        //Log.setLevel(ch.qos.logback.classic.Level.WARN);

	        System.out.println("Starting " + getClass().getSimpleName());
	        this.vmList = new ArrayList<>();
	        this.cloudletList = new ArrayList<>();

	        simulation = new CloudSim();

	        Datacenter datacenter0 = createDatacenter();

	         broker0 = createBroker();

	        createAndSubmitVms(broker0);
	        createAndSubmitCloudlets(broker0);
	        calculateClusterAndFairRatios(cloudletInfo, 1, clusterInfo);
	        Scheduling();
	        simulation.start();

	        List<Cloudlet> finishedCloudlets = broker0.getCloudletFinishedList();
	        new CloudletsTableBuilder(finishedCloudlets).build();
	        printCloudlets();
	        print(broker0);
	    }
		
	   
		//**************************************************************************

		
		private org.cloudbus.cloudsim.brokers.DatacenterBrokerHeuristic createBroker() {
			createSimulatedAnnealingHeuristic();
			org.cloudbus.cloudsim.brokers.DatacenterBrokerHeuristic broker0 = new org.cloudbus.cloudsim.brokers.DatacenterBrokerHeuristic(simulation);
			broker0.setHeuristic(heuristic);
			return broker0;
		}

		private void createAndSubmitCloudlets(org.cloudbus.cloudsim.brokers.DatacenterBrokerHeuristic broker0) {
			int ram_temp =700;
			int cpu_temp =0;
			int k=0;
			for(int j = 0; j < number_Of_Jobs; j++) {
				int cpu_minshare =1; //getRandomNumberOfPes(number_Of_PEs);
				int ram_minshare = 100; //getRandomAmountOfRam(amount_Of_Ram);
				cpu_temp++;
				
			for(int i = 0; i < CLOUDLETS_TO_CREATE; i++){
			   cloudletList.add(createCloudlet(broker0, cpu_temp,ram_temp)) ;
			   cloudletInfo[k][index_cpu_req]= cpu_temp;
			   cloudletInfo[k][index_ram_req]= ram_temp;
			   cloudletInfo[k][index_cpu_minshare]= cpu_minshare;
			   cloudletInfo[k][index_ram_minshare]= ram_minshare;
			   k++;
			  
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
			SA_INITIAL_TEMPERATURE= CalculateTemperature();
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
		//fitness	double roudRobinMappingCost = computeRoudRobinMappingCost();
			  
			printSolution(
			        "Heuristic solution for mapping cloudlets to Vm's         ",
			        heuristic.getBestSolutionSoFar(), false);

			System.out.printf(
			    "The heuristic solution cost represents %.2f%% of the DRF cost used by the DatacenterBrokerSimple\n",
			    heuristic.getBestSolutionSoFar().getCost());
			
			System.out.printf("The solution finding spend %.2f seconds to finish\n", broker0.getHeuristic().getSolveTime());
			System.out.println("Simulated Annealing Parameters");
			System.out.printf("\tInitial Temperature: %.2f", SA_INITIAL_TEMPERATURE);
			System.out.printf(" Cooling Rate: %.4f", SA_COOLING_RATE);
			System.out.printf(" Cold Temperature: %.6f", SA_COLD_TEMPERATURE);
			System.out.printf(" Number of neighborhood searches by iteration: %d\n", SA_NUMBER_OF_NEIGHBORHOOD_SEARCHES);
	        System.out.println(getClass().getSimpleName() + " finished!");
		}
	    //**************************************************************************


		/**
	     * Randomly gets a number of PEs (CPU cores).
	     *
	     * @param maxPesNumber the maximum value to get a random number of PEs
	     * @return the randomly generated PEs number
	     */
	    private int getRandomNumberOfPes(int maxPesNumber) {
	        return heuristic.getRandomValue(maxPesNumber);
	    }
	    //**************************************************************************

	    
	    private int getRandomAmountOfRam(int maxRamAmount) {
	        return heuristic.getRandomValue(maxRamAmount);
	    }
	    //**************************************************************************
  

	    private DatacenterSimple createDatacenter() {
	        List<Host> hostList = new ArrayList<>();
	        for(int i = 0; i < HOSTS_TO_CREATE; i++) {
	            hostList.add(createHost());
	        }

	        return new DatacenterSimple(simulation, hostList, new VmAllocationPolicySimple());
	    }
	    //**************************************************************************


	    private Host createHost() {
	        long mips = 1000; // capacity of each CPU core (in Million Instructions per Second)
	       // int  ram = 8192; // host memory (Megabyte)
	        long storage = 1000000; // host storage
	        long bw = 10000;

	        List<Pe> peList = new ArrayList<>();
	        /*Creates the Host's CPU cores and defines the provisioner
	        used to allocate each core for requesting VMs.*/
	        for(int i = 0; i < 8; i++)
	            peList.add(new PeSimple(mips, new PeProvisionerSimple()));

	       return new HostSimple(amount_Of_Ram, bw, storage, peList)
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
	    								int amountOfRam) {
	        long length = 400000; //in Million Structions (MI)
	        long fileSize = 300; //Size (in bytes) before execution
	        long outputSize = 300; //Size (in bytes) after execution

	        //Defines how CPU, RAM and Bandwidth resources are used
	        //Sets the same utilization model for all these resources.
	        UtilizationModel utilization = new UtilizationModelFull();

	        return new CloudletSimple(numberOfCreatedCloudlets++, length, numberOfPes, amountOfRam)
	        .setFileSize(fileSize)
	        .setOutputSize(outputSize)
	        .setUtilizationModel(utilization);
	    }
	    //**************************************************************************

	    private double Scheduling() {
	             
	           UpdateNeedy((ArrayList<Cloudlet>) cloudletList);
	           
	            if (isNeedy.contains(true))
	            {
	            	//drf
	            	DRfScheduling((ArrayList<Cloudlet>) cloudletList);
	            	
	            }
	            else {
	            	//simulated annealing
	            
	            	SimulatedAnnealingScheduling((ArrayList<Cloudlet>) cloudletList);
	            }
	            
	            
	            
	            
	            
	            //add to usage_cloudlet
	            /*
	            if ( !broker0.getCloudletFinishedList().contains(c.getId())) {
	            usage_Cloudlet[(int) c.getId()][index_cpu_req]+= 
	            		cloudletInfo[(int) c.getId()][index_cpu_req];
	            usage_Cloudlet[(int) c.getId()][index_ram_req]+=
	            		cloudletInfo[(int) c.getId()][index_ram_req];}*/
	        
	        /*
	        printSolution(
	            "Round robin solution used by DatacenterBrokerSimple class",
	            solution, false);
	        return solution.getCost();
	        */
			return 0;
	    }
	  
		//**************************************************************************

	    private double SimulatedAnnealingScheduling(ArrayList<Cloudlet> cloudletList) {
	    	CloudletToVmMappingSolution FairnessSolution =
	                new CloudletToVmMappingSolution(heuristic);
	    	
	    	 Map<Integer, Double> unSortedMap = new HashMap<Integer, Double>() ;
	    	for (Cloudlet c:cloudletList) {
				  //fill sorted but not sorted yet  (cloudlet id, drf shares)
				 
				  
				unSortedMap .put((int)c.getId(), calculateOurFairness(c,1));
			  }
			  
			    //start ascending sort 
			   LinkedHashMap<Integer, Double> sortedMap = new LinkedHashMap<>();
		       unSortedMap.entrySet().stream().sorted(Map.Entry.comparingByValue())
		                .forEachOrdered(x -> sortedMap.put(x.getKey(), x.getValue()));
		       
			  //bind to VM pick with lowest dominant share
		       
		       int i = 0;
		       for (Map.Entry<Integer, Double> entry : sortedMap.entrySet()) {
		    	 //bind the cloudlets to the vms. This way, the broker
		           // will submit the bound cloudlets only to the specific VM
		          
		    	   
		             Cloudlet cloudlet = cloudletList.get(entry.getKey());
		             
		             if (!broker0.getCloudletWaitingList().contains(cloudlet)) {
		            	Vm v=  cloudlet 
		     	                .getBroker()
		     	                .getVmCreatedList()
		     	                .stream()
		     	                .filter(vm -> vm.getNumberOfPes() >= cloudlet.getNumberOfPes())
		     	                .min(Comparator.comparingLong(Vm::getNumberOfPes))
		     	                .orElse(broker0.defaultVmMapper(cloudlet));
		     	                cloudlet.setVm(v);
		     	                AddTOArray(cloudlet);
				        }}
	    	
	    	
	    	
	    	
	       
	        printSolution(
	            "DRF solution used by DatacenterBrokerSimple class",
	            FairnessSolution, false);
	        return FairnessSolution.getCost();
	    }
			
		

		//**************************************************************************

		private void printSolution(String title,
	            CloudletToVmMappingSolution solution,
	            boolean showIndividualCloudletFitness) {
	        System.out.printf("%s (cost %.2f fitness %.6f)\n",
	                title, solution.getCost(), solution.getFitness());
	        if(!showIndividualCloudletFitness)
	            return;

	        for(Map.Entry<Cloudlet, Vm> e: solution.getResult().entrySet()){
	            System.out.printf(
	                "Cloudlet %3d (%d PEs, %6d MI) mapped to Vm %3d (%d PEs, %6.0f MIPS)\n",
	                e.getKey().getId(),
	                e.getKey().getNumberOfPes(), e.getKey().getLength(),
	                e.getValue().getId(),
	                e.getValue().getNumberOfPes(), e.getValue().getMips());
	        }
	        System.out.println();
	    }
		
}
