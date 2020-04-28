clear all
cd ~/Dropbox/MPCwhose


**------------------------------------------------;

**  This program reads the 2008 SIPP Wave 7 Topical Module Data File 
**  Note:  This program is distributed under the GNU GPL. See end of
**  this file and http://www.gnu.org/licenses/ for details.
**  by Jean Roth Fri Feb 17 11:42:26 EST 2012
**  Please report errors to jroth@nber.org
**  run with do sippp08putm7

**-----------------------------------------------;

** The following line should contain
**   the complete path and name of the raw data file.
**   On a PC, use backslashes in paths as in C:\  

global read_dir = "/home/david/Dropbox/MPCwhose/"
** The following line should contain the path to your output '.dta' file 

global save_dir ="/home/david/Dropbox/MPCwhose/"


** The following line should contain the path to the data dictionary file 

local dct_name "sippp08putm7.dct"

** The line below does NOT need to be changed 
local dat_name  "${read_dir}p08putm7.dat"

di "`dat_name'"
quietly infile using "${read_dir}`dct_name'", using(`dat_name') clear


**  Decimal places have been made explict in the dictionary file.
**  Stata resolves a missing value of -1 / # of decimal places as a missing value.

**Everything below this point, aside from the final save, are value labels

#delimit ;

;
label values spanel   spanel; 
label define spanel  
	2008        "Panel Year"                    
;
label values tfipsst  tfipsst;
label define tfipsst 
	1           "Alabama"                       
	2           "Alaska"                        
	4           "Arizona"                       
	5           "Arkansas"                      
	6           "California"                    
	8           "Colorado"                      
	9           "Connecticut"                   
	10          "Delaware"                      
	11          "DC"                            
	12          "Florida"                       
	13          "Georgia"                       
	15          "Hawaii"                        
	16          "Idaho"                         
	17          "Illinois"                      
	18          "Indiana"                       
	19          "Iowa"                          
	20          "Kansas"                        
	21          "Kentucky"                      
	22          "Louisiana"                     
	23          "Maine"                         
	24          "Maryland"                      
	25          "Massachusetts"                 
	26          "Michigan"                      
	27          "Minnesota"                     
	28          "Mississippi"                   
	29          "Missouri"                      
	30          "Montana"                       
	31          "Nebraska"                      
	32          "Nevada"                        
	33          "New Hampshire"                 
	34          "New Jersey"                    
	35          "New Mexico"                    
	36          "New York"                      
	37          "North Carolina"                
	38          "North Dakota"                  
	39          "Ohio"                          
	40          "Oklahoma"                      
	41          "Oregon"                        
	42          "Pennsylvania"                  
	44          "Rhode Island"                  
	45          "South Carolina"                
	46          "South Dakota"                  
	47          "Tennessee"                     
	48          "Texas"                         
	49          "Utah"                          
	50          "Vermont"                       
	51          "Virginia"                      
	53          "Washington"                    
	54          "West Virginia"                 
	55          "Wisconsin"                     
	56          "Wyoming"                       
;
label values eoutcome eoutcome;
label define eoutcome
	201         "Completed interview"           
	203         "Compl. partial- missing data; no"
	207         "Complete partial - TYPE-Z; no" 
	213         "TYPE-A, language problem"      
	216         "TYPE-A, no one home (noh)"     
	217         "TYPE-A, temporarily absent (ta)"
	218         "TYPE-A, hh refused"            
	219         "TYPE-A, other occupied (specify)"
	234         "TYPE-B, entire hh institut. or"
	248         "TYPE-C, other (specify)"       
	249         "TYPE-C, sample adjustment"     
	250         "TYPE-C, hh deceased"           
	251         "TYPE-C, moved out of country"  
	252         "TYPE-C, living in armed forces"
	253         "TYPE-C, on active duty in Armed"
	254         "TYPE-C, no one over age 15 years"
	255         "TYPE-C, no Wave 1 persons"     
	260         "TYPE-D, moved address unknown" 
	261         "TYPE-D, moved within U.S. but" 
	262         "TYPE-C, merged with another SIPP"
	270         "TYPE-C, mover, no longer located"
	271         "TYPE-C, mover, new address"    
	280         "TYPE-D, mover, no longer located"
;
label values rfid2    rfid2l; 
label define rfid2l  
	-1          "Not in Universe"               
;
label values epopstat epopstat;
label define epopstat
	1           "Adult (15 years of age or older)"
	2           "Child (Under 15 years of age)" 
;
label values eppintvw eppintvw;
label define eppintvw
	1           "Interview (self)"              
	2           "Interview (proxy)"             
	3           "Noninterview - Type Z"         
	4           "Noninterview - pseudo Type Z." 
	5           "Children under 15 during"      
;
label values esex     esex;   
label define esex    
	1           "Male"                          
	2           "Female"                        
;
label values erace    erace;  
label define erace   
	1           "White alone"                   
	2           "Black alone"                   
	3           "Asian alone"                   
	4           "Residual"                      
;
label values eorigin  eorigin;
label define eorigin 
	1           "Yes"                           
	2           "No"                            
;
label values errp     errp;   
label define errp    
	1           "Reference person with related" 
	2           "Reference Person without related"
	3           "Spouse of reference person"    
	4           "Child of reference person"     
	5           "Grandchild of reference person"
	6           "Parent of reference person"    
	7           "Brother/sister of reference person"
	8           "Other relative of reference person"
	9           "Foster child of reference person"
	10          "Unmarried partner of reference"
	11          "Housemate/roommate"            
	12          "Roomer/boarder"                
	13          "Other non-relative of reference"
;
label values tage     tage;   
label define tage    
	0           "Less than 1 full year old"     
;
label values ems      ems;    
label define ems     
	1           "Married, spouse present"       
	2           "Married, spouse absent"        
	3           "Widowed"                       
	4           "Divorced"                      
	5           "Separated"                     
	6           "Never Married"                 
;
label values epnspous epnspous;
label define epnspous
	9999        "Spouse not in household or person"
;
label values epnmom   epnmom; 
label define epnmom  
	9999        "No mother in household"        
;
label values epndad   epndad; 
label define epndad  
	9999        "No father in household"        
;
label values epnguard epnguard;
label define epnguard
	-1          "Not in Universe"               
	9999        "Guardian not in household"     
;
label values rdesgpnt rdesgpnt;
label define rdesgpnt
	-1          "Not in Universe"               
	1           "Yes"                           
	2           "No"                            
;
label values eeducate eeducate;
label define eeducate
	-1          "Not in Universe"               
	31          "Less Than 1st Grade"           
	32          "1st, 2nd, 3rd or 4th grade"    
	33          "5th Or 6th Grade"              
	34          "7th Or 8th Grade"              
	35          "9th Grade"                     
	36          "10th Grade"                    
	37          "11th Grade"                    
	38          "12th grade, no diploma"        
	39          "High School Graduate - (diploma"
	40          "Some college, but no degree"   
	41          "Diploma or certificate from a" 
	43          "Associate (2-yr) college degree"
	44          "Bachelor's degree (for example:"
	45          "Master's degree (For example: MA,"
	46          "Professional School degree (for"
	47          "Doctorate degree (for example:"
;
label values sinthhid sinthhid;
label define sinthhid
	0           "Not In Universe"               
;
label values ealunv   ealunv; 
label define ealunv  
	-1          "Not in Universe"               
	1           "In universe"                   
;
label values ealr     ealr;   
label define ealr    
	-1          "Not in Universe"               
	1           "Yes"                           
	2           "No"                            
;
label values aalr     aalr;   
label define aalr    
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values ealry    ealry;  
label define ealry   
	-1          "Not in Universe"               
;
label values aalry    aalry;  
label define aalry   
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values talrb    talrb;  
label define talrb   
	0           "None or not in universe"       
;
label values aalrb    aalrb;  
label define aalrb   
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values ealra1   ealra1l;
label define ealra1l 
	-1          "Not in Universe"               
	1           "Certificates of deposit or other"
	2           "Money market funds"            
	3           "U.S. Government securities"    
	4           "Municipal or corporate bonds"  
	5           "U.S. Savings Bonds"            
	6           "Stocks or mutual fund shares"  
	7           "Other assets"                  
;
label values aalra1   aalra1l;
label define aalra1l 
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values ealra2   ealra2l;
label define ealra2l 
	-1          "Not in Universe"               
	1           "Certificates of deposit or other"
	2           "Money market funds"            
	3           "U.S. Government securities"    
	4           "Municipal or corporate bonds"  
	5           "U.S. Savings Bonds"            
	6           "Stocks or mutual fund shares"  
	7           "Other assets"                  
;
label values aalra2   aalra2l;
label define aalra2l 
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values ealra3   ealra3l;
label define ealra3l 
	-1          "Not in Universe"               
	1           "Certificates of deposit or other"
	2           "Money market funds"            
	3           "U.S. Government securities"    
	4           "Municipal or corporate bonds"  
	5           "U.S. Savings Bonds"            
	6           "Stocks or mutual fund shares"  
	7           "Other assets"                  
;
label values aalra3   aalra3l;
label define aalra3l 
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values ealra4   ealra4l;
label define ealra4l 
	-1          "Not in Universe"               
	1           "Certificates of deposit or other"
	2           "Money market funds"            
	3           "U.S. Government securities"    
	4           "Municipal or corporate bonds"  
	5           "U.S. Savings Bonds"            
	6           "Stocks or mutual fund shares"  
	7           "Other assets"                  
;
label values aalra4   aalra4l;
label define aalra4l 
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values ealk     ealk;   
label define ealk    
	-1          "Not in Universe"               
	1           "Yes"                           
	2           "No"                            
;
label values aalk     aalk;   
label define aalk    
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values ealky    ealky;  
label define ealky   
	-1          "Not in Universe"               
;
label values aalky    aalky;  
label define aalky   
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values talkb    talkb;  
label define talkb   
	0           "None or not in universe"       
;
label values aalkb    aalkb;  
label define aalkb   
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values ealka1   ealka1l;
label define ealka1l 
	-1          "Not in Universe"               
	1           "Certificates of deposit or other"
	2           "Money market funds"            
	3           "U.S. Government securities"    
	4           "Municipal or corporate bonds"  
	5           "U.S. Savings Bonds"            
	6           "Stocks or mutual fund shares"  
	7           "Other assets"                  
;
label values aalka1   aalka1l;
label define aalka1l 
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values ealka2   ealka2l;
label define ealka2l 
	-1          "Not in Universe"               
	1           "Certificates of deposit or other"
	2           "Money market funds"            
	3           "U.S. Government securities"    
	4           "Municipal or corporate bonds"  
	5           "U.S. Savings Bonds"            
	6           "Stocks or mutual fund shares"  
	7           "Other assets"                  
;
label values aalka2   aalka2l;
label define aalka2l 
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values ealka3   ealka3l;
label define ealka3l 
	-1          "Not in Universe"               
	1           "Certificates of deposit or other"
	2           "Money market funds"            
	3           "U.S. Government securities"    
	4           "Municipal or corporate bonds"  
	5           "U.S. Savings Bonds"            
	6           "Stocks or mutual fund shares"  
	7           "Other assets"                  
;
label values aalka3   aalka3l;
label define aalka3l 
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values ealka4   ealka4l;
label define ealka4l 
	-1          "Not in Universe"               
	1           "Certificates of deposit or other"
	2           "Money market funds"            
	3           "U.S. Government securities"    
	4           "Municipal or corporate bonds"  
	5           "U.S. Savings Bonds"            
	6           "Stocks or mutual fund shares"  
	7           "Other assets"                  
;
label values aalka4   aalka4l;
label define aalka4l 
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values ealt     ealt;   
label define ealt    
	-1          "Not in Universe"               
	1           "Yes"                           
	2           "No"                            
;
label values aalt     aalt;   
label define aalt    
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values ealty    ealty;  
label define ealty   
	-1          "Not in Universe"               
;
label values aalty    aalty;  
label define aalty   
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values taltb    taltb;  
label define taltb   
	0           "None or not in universe"       
;
label values aaltb    aaltb;  
label define aaltb   
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values ealta1   ealta1l;
label define ealta1l 
	-1          "Not in Universe"               
	1           "Certificates of deposit or other"
	2           "Money market funds"            
	3           "U.S. Government securities"    
	4           "Municipal or corporate bonds"  
	5           "U.S. Savings Bonds"            
	6           "Stocks or mutual fund shares"  
	7           "Other assets"                  
;
label values aalta1   aalta1l;
label define aalta1l 
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values ealta2   ealta2l;
label define ealta2l 
	-1          "Not in Universe"               
	1           "Certificates of deposit or other"
	2           "Money market funds"            
	3           "U.S. Government securities"    
	4           "Municipal or corporate bonds"  
	5           "U.S. Savings Bonds"            
	6           "Stocks or mutual fund shares"  
	7           "Other assets"                  
;
label values aalta2   aalta2l;
label define aalta2l 
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values ealta3   ealta3l;
label define ealta3l 
	-1          "Not in Universe"               
	1           "Certificates of deposit or other"
	2           "Money market funds"            
	3           "U.S. Government securities"    
	4           "Municipal or corporate bonds"  
	5           "U.S. Savings Bonds"            
	6           "Stocks or mutual fund shares"  
	7           "Other assets"                  
;
label values aalta3   aalta3l;
label define aalta3l 
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values ealta4   ealta4l;
label define ealta4l 
	-1          "Not in Universe"               
	1           "Certificates of deposit or other"
	2           "Money market funds"            
	3           "U.S. Government securities"    
	4           "Municipal or corporate bonds"  
	5           "U.S. Savings Bonds"            
	6           "Stocks or mutual fund shares"  
	7           "Other assets"                  
;
label values aalta4   aalta4l;
label define aalta4l 
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values ealow    ealow;  
label define ealow   
	-1          "Not in Universe"               
	1           "Yes"                           
	2           "No"                            
;
label values aalow    aalow;  
label define aalow   
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values talowa   talowa; 
label define talowa  
	0           "Not In Universe"               
;
label values aalowa   aalowa; 
label define aalowa  
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values ealsb    ealsb;  
label define ealsb   
	-1          "Not in Universe"               
	1           "Yes"                           
	2           "No"                            
;
label values aalsb    aalsb;  
label define aalsb   
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values talsbv   talsbv; 
label define talsbv  
	0           "Not In Universe"               
;
label values aalsbv   aalsbv; 
label define aalsbv  
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values ealjch   ealjch; 
label define ealjch  
	-1          "Not in Universe"               
	1           "Yes"                           
	2           "No"                            
;
label values aaljch   aaljch; 
label define aaljch  
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values taljcha  taljcha;
label define taljcha 
	0           "None or not in universe"       
;
label values aaljcha  aaljcha;
label define aaljcha 
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values ealjdb   ealjdb; 
label define ealjdb  
	-1          "Not in Universe"               
	1           "Yes"                           
	2           "No"                            
;
label values aaljdb   aaljdb; 
label define aaljdb  
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values ealjdl   ealjdl; 
label define ealjdl  
	-1          "Not in Universe"               
	1           "Yes"                           
	2           "No"                            
;
label values aaljdl   aaljdl; 
label define aaljdl  
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values ealjdo   ealjdo; 
label define ealjdo  
	-1          "Not in Universe"               
	1           "Yes"                           
	2           "No"                            
;
label values aaljdo   aaljdo; 
label define aaljdo  
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values taljdab  taljdab;
label define taljdab 
	0           "Not In Universe"               
;
label values aaljdab  aaljdab;
label define aaljdab 
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values taljdal  taljdal;
label define taljdal 
	0           "Not In Universe"               
;
label values aaljdal  aaljdal;
label define aaljdal 
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values taljdao  taljdao;
label define taljdao 
	0           "Not In Universe"               
;
label values aaljdao  aaljdao;
label define aaljdao 
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values ealich   ealich; 
label define ealich  
	-1          "Not in Universe"               
	1           "Yes"                           
	2           "No"                            
;
label values aalich   aalich; 
label define aalich  
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values talicha  talicha;
label define talicha 
	0           "None or not in universe"       
;
label values aalicha  aalicha;
label define aalicha 
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values ealil    ealil;  
label define ealil   
	-1          "Not in Universe"               
	1           "Yes"                           
	2           "No"                            
;
label values aalil    aalil;  
label define aalil   
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values ealidb   ealidb; 
label define ealidb  
	-1          "Not in Universe"               
	1           "Yes"                           
	2           "No"                            
;
label values aalidb   aalidb; 
label define aalidb  
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values ealidl   ealidl; 
label define ealidl  
	-1          "Not in Universe"               
	1           "Yes"                           
	2           "No"                            
;
label values aalidl   aalidl; 
label define aalidl  
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values ealido   ealido; 
label define ealido  
	-1          "Not in Universe"               
	1           "Yes"                           
	2           "No"                            
;
label values aalido   aalido; 
label define aalido  
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values talidab  talidab;
label define talidab 
	0           "Not In Universe"               
;
label values aalidab  aalidab;
label define aalidab 
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values talidal  talidal;
label define talidal 
	0           "Not In Universe"               
;
label values aalidal  aalidal;
label define aalidal 
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values talidao  talidao;
label define talidao 
	0           "Not In Universe"               
;
label values aalidao  aalidao;
label define aalidao 
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values ealli    ealli;  
label define ealli   
	-1          "Not in Universe"               
	1           "Yes"                           
	2           "No"                            
;
label values aalli    aalli;  
label define aalli   
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values talliv   talliv; 
label define talliv  
	0           "Zero or not in universe"       
;
label values aalliv   aalliv; 
label define aalliv  
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values eallit   eallit; 
label define eallit  
	-1          "Not in Universe"               
	1           "Term only"                     
	2           "Whole life only"               
	3           "Both types"                    
;
label values aallit   aallit; 
label define aallit  
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values eallie   eallie; 
label define eallie  
	-1          "Not in Universe"               
	1           "Yes"                           
	2           "No"                            
;
label values aallie   aallie; 
label define aallie  
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values talliev  talliev;
label define talliev 
	0           "Zero or not in universe"       
;
label values aalliev  aalliev;
label define aalliev 
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values ehreunv  ehreunv;
label define ehreunv 
	-1          "Not in Universe"               
	1           "In universe"                   
;
label values eremobho eremobho;
label define eremobho
	-1          "Not in Universe"               
	1           "Yes"                           
	2           "No"                            
;
label values aremobho aremobho;
label define aremobho
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values ehowner1 ehownero;
label define ehownero
	-1          "Not in Universe"               
;
label values ahowner1 ahownero;
label define ahownero
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values ehowner2 ehownerk;
label define ehownerk
	-1          "Not in Universe"               
;
label values ahowner2 ahownerk;
label define ahownerk
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)}"
;
label values ehowner3 ehownerl;
label define ehownerl
	-1          "Not in Universe"               
;
label values ehbuymo  ehbuymo;
label define ehbuymo 
	-1          "Not in Universe"               
;
label values ahbuymo  ahbuymo;
label define ahbuymo 
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values ehbuyyr  ehbuyyr;
label define ehbuyyr 
	-1          "Not in Universe"               
;
label values ahbuyyr  ahbuyyr;
label define ahbuyyr 
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values ehmort   ehmort; 
label define ehmort  
	-1          "Not in Universe"               
	1           "Yes"                           
	2           "No"                            
;
label values ahmort   ahmort; 
label define ahmort  
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values enummort enummort;
label define enummort
	-1          "Not in Universe"               
;
label values anummort anummort;
label define anummort
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values tmor1pr  tmor1pr;
label define tmor1pr 
	0           "Not In Universe"               
;
label values amor1pr  amor1pr;
label define amor1pr 
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values emor1yr  emor1yr;
label define emor1yr 
	-1          "Not in Universe"               
;
label values amor1yr  amor1yr;
label define amor1yr 
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values emor1mo  emor1mo;
label define emor1mo 
	-1          "Not in Universe"               
;
label values amor1mo  amor1mo;
label define amor1mo 
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values tmor1amt tmor1amt;
label define tmor1amt
	0           "None or not in universe"       
;
label values amor1amt amor1amt;
label define amor1amt
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values tmor1yrs tmor1yrs;
label define tmor1yrs
	-1          "Not in Universe"               
;
label values amor1yrs amor1yrs;
label define amor1yrs
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values emor1int emor1int;
label define emor1int
	-1          "Not in Universe"               
;
label values amor1int amor1int;
label define amor1int
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values emor1var emor1var;
label define emor1var
	-1          "Not in Universe"               
	1           "Variable interest rate"        
	2           "Fixed interest rate"           
;
label values amor1var amor1var;
label define amor1var
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values emor1pgm emor1pgm;
label define emor1pgm
	-1          "Not in Universe"               
	1           "Yes - FHA LOAN"                
	2           "Yes - VA LOAN"                 
	3           "NO"                            
;
label values amor1pgm amor1pgm;
label define amor1pgm
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values tmor2pr  tmor2pr;
label define tmor2pr 
	0           "None or not in universe"       
	1           "Flag indicating principal on"  
;
label values amor2pr  amor2pr;
label define amor2pr 
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values emor2yr  emor2yr;
label define emor2yr 
	-1          "Not in Universe"               
;
label values amor2yr  amor2yr;
label define amor2yr 
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values emor2mo  emor2mo;
label define emor2mo 
	-1          "Not in Universe"               
;
label values amor2mo  amor2mo;
label define amor2mo 
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values tmor2amt tmor2amt;
label define tmor2amt
	0           "None or not in universe"       
	1           "Flag indicating reported amount"
;
label values amor2amt amor2amt;
label define amor2amt
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values tmor2yrs tmor2yrs;
label define tmor2yrs
	-1          "Not in Universe"               
;
label values amor2yrs amor2yrs;
label define amor2yrs
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values emor2int emor2int;
label define emor2int
	-1          "Not in Universe"               
;
label values amor2int amor2int;
label define amor2int
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values emor2var emor2var;
label define emor2var
	-1          "Not in Universe"               
	1           "Variable interest rate"        
	2           "Fixed interest rate"           
;
label values amor2var amor2var;
label define amor2var
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values emor2pgm emor2pgm;
label define emor2pgm
	-1          "Not in Universe"               
	1           "Yes-FHA LOAN"                  
	2           "Yes-VA LOAN"                   
	3           "NO"                            
;
label values amor2pgm amor2pgm;
label define amor2pgm
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values tmor3pr  tmor3pr;
label define tmor3pr 
	0           "None or not in universe"       
	1           "Flag indicating principal reported"
;
label values amor3pr  amor3pr;
label define amor3pr 
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values tpropval tpropval;
label define tpropval
	0           "None or not in universe"       
;
label values apropval apropval;
label define apropval
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values emhloan  emhloan;
label define emhloan 
	-1          "Not in Universe"               
	1           "Yes"                           
	2           "No"                            
;
label values amhloan  amhloan;
label define amhloan 
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values emhtype  emhtype;
label define emhtype 
	-1          "Not in Universe"               
	1           "Mobile home only"              
	2           "Site only"                     
	3           "Site and home"                 
;
label values amhtype  amhtype;
label define amhtype 
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values tmhpr    tmhpr;  
label define tmhpr   
	0           "None or not in universe"       
;
label values amhpr    amhpr;  
label define amhpr   
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values tmhval   tmhval; 
label define tmhval  
	0           "None or not in universe"       
;
label values amhval   amhval; 
label define amhval  
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values thomeamt thomeamt;
label define thomeamt
	0           "None or not in universe"       
;
label values ahomeamt ahomeamt;
label define ahomeamt
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values tutils   tutils; 
label define tutils  
	0           "None or not in universe"       
;
label values autils   autils; 
label define autils  
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values eperspay eperspay;
label define eperspay
	-1          "Not in Universe"               
	1           "Yes"                           
	2           "No"                            
;
label values aperspay aperspay;
label define aperspay
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values eperspya eperspya;
label define eperspya
	-1          "Not in Universe"               
;
label values aperspya aperspya;
label define aperspya
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values eperspy1 eperspyk;
label define eperspyk
	-1          "Not in Universe"               
;
label values aperspy1 aperspyk;
label define aperspyk
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values eperspy2 eperspyl;
label define eperspyl
	-1          "Not in Universe"               
;
label values eperspy3 eperspym;
label define eperspym
	-1          "Not in Universe"               
;
label values tpersam1 tpersama;
label define tpersama
	0           "None or not in universe"       
;
label values apersam1 apersama;
label define apersama
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values tpersam2 tpersamk;
label define tpersamk
	0           "None or not in universe"       
;
label values apersam2 apersamk;
label define apersamk
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values tpersam3 tpersaml;
label define tpersaml
	0           "None or not in universe"       
;
label values apersam3 apersaml;
label define apersaml
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values epaycare epaycare;
label define epaycare
	-1          "Not in Universe"               
	1           "Yes"                           
	2           "No"                            
;
label values apaycare apaycare;
label define apaycare
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values tcarecst tcarecst;
label define tcarecst
	0           "None or not in universe"       
;
label values acarecst acarecst;
label define acarecst
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values eothre   eothre; 
label define eothre  
	-1          "Not in Universe"               
	1           "Yes"                           
	2           "No"                            
;
label values aothre   aothre; 
label define aothre  
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values eothreo1 eothreoe;
label define eothreoe
	-1          "Not in Universe"               
;
label values aothreo1 aothreoe;
label define aothreoe
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values eothreo2 eothreok;
label define eothreok
	-1          "Not in Universe"               
;
label values eothreo3 eothreol;
label define eothreol
	-1          "Not in Universe"               
;
label values tothreva tothreva;
label define tothreva
	0           "None or not in universe"       
;
label values aothreva aothreva;
label define aothreva
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values eautoown eautoown;
label define eautoown
	-1          "Not in Universe"               
	1           "Yes"                           
	2           "No"                            
;
label values aautoown aautoown;
label define aautoown
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values eautonum eautonum;
label define eautonum
	-1          "Not in Universe"               
;
label values aautonum aautonum;
label define aautonum
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values ea1own1  ea1own1l;
label define ea1own1l
	-1          "Not in Universe"               
;
label values aa1own1  aa1own1l;
label define aa1own1l
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values ea1own2  ea1own2l;
label define ea1own2l
	-1          "Not in Universe"               
;
label values tcarval1 tcarvalm;
label define tcarvalm
	0           "None or not in universe"       
;
label values acarval1 acarvalm;
label define acarvalm
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values ta1year  ta1year;
label define ta1year 
	-1          "Not in Universe"               
	9999        "Don't Know, Refusal, Blanks from"
;
label values ea1owed  ea1owed;
label define ea1owed 
	-1          "Not in Universe"               
	1           "Money owed"                    
	2           "Free and clear"                
;
label values aa1owed  aa1owed;
label define aa1owed 
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values ta1amt   ta1amt; 
label define ta1amt  
	0           "None or not in universe"       
;
label values aa1amt   aa1amt; 
label define aa1amt  
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values ea1use   ea1use; 
label define ea1use  
	-1          "Not in Universe"               
	1           "Yes"                           
	2           "No"                            
;
label values aa1use   aa1use; 
label define aa1use  
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values ea2own1  ea2own1l;
label define ea2own1l
	-1          "Not in Universe"               
;
label values aa2own1  aa2own1l;
label define aa2own1l
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values ea2own2  ea2own2l;
label define ea2own2l
	-1          "Not in Universe"               
;
label values tcarval2 tcarvale;
label define tcarvale
	0           "None or not in universe"       
;
label values acarval2 acarvale;
label define acarvale
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values ta2year  ta2year;
label define ta2year 
	-1          "Not in Universe"               
	1986        "Recode for year less than 1986"
	1991        "Recode for year 1987-1991"     
	9999        "Don't Know, Refusal, Blanks from"
;
label values ea2owed  ea2owed;
label define ea2owed 
	-1          "Not in Universe"               
	1           "Money owed"                    
	2           "Free and clear"                
;
label values aa2owed  aa2owed;
label define aa2owed 
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values ta2amt   ta2amt; 
label define ta2amt  
	0           "None or not in universe"       
;
label values aa2amt   aa2amt; 
label define aa2amt  
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values ea2use   ea2use; 
label define ea2use  
	-1          "Not in Universe"               
	1           "Yes"                           
	2           "No"                            
;
label values aa2use   aa2use; 
label define aa2use  
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values ea3own1  ea3own1l;
label define ea3own1l
	-1          "Not in Universe"               
;
label values aa3own1  aa3own1l;
label define aa3own1l
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values ea3own2  ea3own2l;
label define ea3own2l
	-1          "Not in Universe"               
;
label values tcarval3 tcarvalk;
label define tcarvalk
	0           "None or not in universe"       
;
label values acarval3 acarvalk;
label define acarvalk
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values ta3year  ta3year;
label define ta3year 
	-1          "Not in Universe"               
	1971        "Recode for year less than 1971"
	1980        "Recode for year 1972-1980"     
	1986        "Recode for year 1981-1986"     
	1989        "Recode for year 1987-1989"     
	1991        "Recode for year 1990-1991"     
	9999        "Don't Know, Refusal, Blanks from"
;
label values ea3owed  ea3owed;
label define ea3owed 
	-1          "Not in Universe"               
	1           "Money owed"                    
	2           "Free and clear"                
;
label values aa3owed  aa3owed;
label define aa3owed 
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values ta3amt   ta3amt; 
label define ta3amt  
	0           "None or not in universe"       
;
label values aa3amt   aa3amt; 
label define aa3amt  
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values ea3use   ea3use; 
label define ea3use  
	-1          "Not in Universe"               
	1           "Yes"                           
	2           "No"                            
;
label values aa3use   aa3use; 
label define aa3use  
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values eothveh  eothveh;
label define eothveh 
	-1          "Not in Universe"               
	1           "Yes"                           
	2           "No"                            
;
label values aothveh  aothveh;
label define aothveh 
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values eovmtrcy eovmtrcy;
label define eovmtrcy
	-1          "Not in Universe"               
	1           "Yes"                           
	2           "No"                            
;
label values aovmtrcy aovmtrcy;
label define aovmtrcy
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values eovboat  eovboat;
label define eovboat 
	-1          "Not in Universe"               
	1           "Yes"                           
	2           "No"                            
;
label values aovboat  aovboat;
label define aovboat 
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values eovrv    eovrv;  
label define eovrv   
	-1          "Not in Universe"               
	1           "Yes"                           
	2           "Not"                           
;
label values aovrv    aovrv;  
label define aovrv   
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values eovothrv eovothrv;
label define eovothrv
	-1          "Not in Universe"               
	1           "Yes"                           
	2           "Not"                           
;
label values aovothrv aovothrv;
label define aovothrv
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values eov1own1 eov1ownv;
label define eov1ownv
	-1          "Not in Universe"               
;
label values aov1own1 aov1ownv;
label define aov1ownv
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values eov1own2 eov1ownk;
label define eov1ownk
	-1          "Not in Universe"               
;
label values tov1val  tov1val;
label define tov1val 
	0           "None or not in universe"       
;
label values aov1val  aov1val;
label define aov1val 
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values eov1owe  eov1owe;
label define eov1owe 
	-1          "Not in Universe"               
	1           "Money owed"                    
	2           "Free and clear"                
;
label values aov1owe  aov1owe;
label define aov1owe 
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values tov1amt  tov1amt;
label define tov1amt 
	0           "None or not in universe"       
;
label values aov1amt  aov1amt;
label define aov1amt 
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values eov2own1 eov2ownt;
label define eov2ownt
	-1          "Not in Universe"               
;
label values aov2own1 aov2ownt;
label define aov2ownt
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values eov2own2 eov2ownk;
label define eov2ownk
	-1          "Not in Universe"               
;
label values tov2val  tov2val;
label define tov2val 
	0           "None or not in universe"       
;
label values aov2val  aov2val;
label define aov2val 
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values eov2owe  eov2owe;
label define eov2owe 
	-1          "Not in Universe"               
	1           "Money owed"                    
	2           "Free and clear"                
;
label values aov2owe  aov2owe;
label define aov2owe 
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values tov2amt  tov2amt;
label define tov2amt 
	0           "None or not in universe"       
;
label values aov2amt  aov2amt;
label define aov2amt 
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values thhtnw   thhtnw; 
label define thhtnw  
	0           "None or not in universe"       
;
label values thhtwlth thhtwlth;
label define thhtwlth
	0           "None or not in universe"       
;
label values thhtheq  thhtheq;
label define thhtheq 
	0           "None or not in universe"       
;
label values thhmortg thhmortg;
label define thhmortg
	0           "None or not in universe"       
;
label values thhvehcl thhvehcl;
label define thhvehcl
	0           "None or not in universe"       
;
label values thhbeq   thhbeq; 
label define thhbeq  
	0           "None or not in universe"       
;
label values thhintbk thhintbk;
label define thhintbk
	0           "None or not in universe"       
;
label values thhintot thhintot;
label define thhintot
	0           "None or not in universe"       
;
label values thhstk   thhstk; 
label define thhstk  
	0           "None or not in universe"       
;
label values thhore   thhore; 
label define thhore  
	0           "None or not in universe"       
;
label values thhotast thhotast;
label define thhotast
	0           "None or not in universe"       
;
label values thhira   thhira; 
label define thhira  
	0           "None or not in universe"       
;
label values thhthrif thhthrif;
label define thhthrif
	0           "None or not in universe"       
;
label values thhdebt  thhdebt;
label define thhdebt 
	0           "None or not in universe"       
;
label values thhscdbt thhscdbt;
label define thhscdbt
	0           "None or not in universe"       
;
label values thhuscbt thhuscbt;
label define thhuscbt
	0           "None or not in universe"       
;
label values eaoaunv  eaoaunv;
label define eaoaunv 
	-1          "Not in Universe"               
	1           "In universe"                   
;
label values toaeq    toaeq;  
label define toaeq   
	0           "None or not in universe"       
;
label values aoaeq    aoaeq;  
label define aoaeq   
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values tiajta   tiajta; 
label define tiajta  
	0           "None or not in universe"       
;
label values aiajta   aiajta; 
label define aiajta  
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values tiaita   tiaita; 
label define tiaita  
	0           "None or not in universe"       
;
label values aiaita   aiaita; 
label define aiaita  
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values timja    timja;  
label define timja   
	0           "None or not in universe"       
;
label values aimja    aimja;  
label define aimja   
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values timia    timia;  
label define timia   
	0           "None or not in universe"       
;
label values aimia    aimia;  
label define aimia   
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values esmjm    esmjm;  
label define esmjm   
	-1          "Not in Universe"               
	1           "Yes"                           
	2           "No"                            
;
label values asmjm    asmjm;  
label define asmjm   
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values esmjs    esmjs;  
label define esmjs   
	-1          "Not in Universe"               
	1           "Yes"                           
	2           "No"                            
;
label values asmjs    asmjs;  
label define asmjs   
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values tsmjv    tsmjv;  
label define tsmjv   
	0           "None or not in universe"       
;
label values asmjv    asmjv;  
label define asmjv   
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values esmjma   esmjma; 
label define esmjma  
	-1          "Not in Universe"               
	1           "Yes"                           
	2           "No"                            
;
label values asmjma   asmjma; 
label define asmjma  
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values tsmjmav  tsmjmav;
label define tsmjmav 
	0           "None or not in universe"       
;
label values asmjmav  asmjmav;
label define asmjmav 
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values esmi     esmi;   
label define esmi    
	-1          "Not in Universe"               
	1           "Yes"                           
	2           "No"                            
;
label values asmi     asmi;   
label define asmi    
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values tsmiv    tsmiv;  
label define tsmiv   
	0           "None or not in universe"       
;
label values asmiv    asmiv;  
label define asmiv   
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values esmima   esmima; 
label define esmima  
	-1          "Not in Universe"               
	1           "Yes"                           
	2           "No"                            
;
label values asmima   asmima; 
label define asmima  
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values tsmimav  tsmimav;
label define tsmimav 
	0           "None or not in universe"       
;
label values asmimav  asmimav;
label define asmimav 
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values erjown   erjown; 
label define erjown  
	-1          "Not in Universe"               
	1           "Yes"                           
	2           "No"                            
;
label values arjown   arjown; 
label define arjown  
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values erjnum   erjnum; 
label define erjnum  
	0           "None or not in universe"       
;
label values arjnum   arjnum; 
label define arjnum  
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values erjtyp1  erjtyp1l;
label define erjtyp1l
	-1          "Not in Universe"               
	1           "Vacation home"                 
	2           "Other residential property"    
	3           "Farm property"                 
	4           "Commercial property"           
	5           "Equipment"                     
	6           "Other"                         
;
label values arjtyp1  arjtyp1l;
label define arjtyp1l
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values erjtyp2  erjtyp2l;
label define erjtyp2l
	-1          "Not in Universe"               
	1           "Vacation home"                 
	2           "Other residential property"    
	3           "Farm property"                 
	4           "Commercial property"           
	5           "Equipment"                     
	6           "Other"                         
;
label values arjtyp2  arjtyp2l;
label define arjtyp2l
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values erjtyp3  erjtyp3l;
label define erjtyp3l
	-1          "Not in Universe"               
	1           "Vacation home"                 
	2           "Other residential property"    
	3           "Farm property"                 
	4           "Commercial property"           
	5           "Equipment"                     
	6           "Other"                         
;
label values arjtyp3  arjtyp3l;
label define arjtyp3l
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values erjtyp4  erjtyp4l;
label define erjtyp4l
	-1          "Not in Universe"               
	1           "Vacation home"                 
	2           "Other residential property"    
	3           "Farm property"                 
	4           "Commercial property"           
	5           "Equipment"                     
	6           "Other"                         
;
label values arjtyp4  arjtyp4l;
label define arjtyp4l
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values erjtyp5  erjtyp5l;
label define erjtyp5l
	-1          "Not in Universe"               
	1           "Vacation home"                 
	2           "Other residential property"    
	3           "Farm property"                 
	4           "Commercial property"           
	5           "Equipment"                     
	6           "Other"                         
;
label values arjtyp5  arjtyp5l;
label define arjtyp5l
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values erjtyp6  erjtyp6l;
label define erjtyp6l
	-1          "Not in Universe"               
	1           "Vacation home"                 
	2           "Other residential property"    
	3           "Farm property"                 
	4           "Commercial property"           
	5           "Equipment"                     
	6           "Other"                         
;
label values arjtyp6  arjtyp6l;
label define arjtyp6l
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values erjat    erjat;  
label define erjat   
	-1          "Not in Universe"               
	1           "Yes"                           
	2           "No"                            
;
label values arjat    arjat;  
label define arjat   
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values erjata   erjata; 
label define erjata  
	-1          "Not in Universe"               
	1           "Yes"                           
	2           "No"                            
;
label values arjata   arjata; 
label define arjata  
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values trjmv    trjmv;  
label define trjmv   
	0           "None or not in universe"       
;
label values arjmv    arjmv;  
label define arjmv   
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values erjdeb   erjdeb; 
label define erjdeb  
	-1          "Not in Universe"               
	1           "Yes"                           
	2           "No"                            
;
label values arjdeb   arjdeb; 
label define arjdeb  
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values trjpri   trjpri; 
label define trjpri  
	0           "None or not in universe"       
;
label values arjpri   arjpri; 
label define arjpri  
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values eriown   eriown; 
label define eriown  
	-1          "Not in Universe"               
	1           "Yes"                           
	2           "No"                            
;
label values ariown   ariown; 
label define ariown  
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values erinum   erinum; 
label define erinum  
	0           "None or not in universe"       
;
label values arinum   arinum; 
label define arinum  
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values eritype1 eritypem;
label define eritypem
	-1          "Not in Universe"               
	1           "Vacation home"                 
	2           "Other residential property"    
	3           "Farm property"                 
	4           "Commercial property"           
	5           "Equipment"                     
	6           "Other"                         
;
label values aritype1 aritypem;
label define aritypem
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values eritype2 eritypek;
label define eritypek
	-1          "Not in Universe"               
	1           "Vacation home"                 
	2           "Other residential property"    
	3           "Farm property"                 
	4           "Commercial property"           
	5           "Equipment"                     
	6           "Other"                         
;
label values aritype2 aritypek;
label define aritypek
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values eritype3 eritypel;
label define eritypel
	-1          "Not in Universe"               
	1           "Vacation home"                 
	2           "Other residential property"    
	3           "Farm property"                 
	4           "Commercial property"           
	5           "Equipment"                     
	6           "Other"                         
;
label values aritype3 aritypel;
label define aritypel
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values eritype4 eritypen;
label define eritypen
	-1          "Not in Universe"               
	1           "Vacation home"                 
	2           "Other residential property"    
	3           "Farm property"                 
	4           "Commercial property"           
	5           "Equipment"                     
	6           "Other"                         
;
label values aritype4 aritypen;
label define aritypen
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values eritype5 eritypeo;
label define eritypeo
	-1          "Not in Universe"               
	1           "Vacation home"                 
	2           "Other residential property"    
	3           "Farm property"                 
	4           "Commercial property"           
	5           "Equipment"                     
	6           "Other"                         
;
label values aritype5 aritypeo;
label define aritypeo
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values eritype6 eritypep;
label define eritypep
	-1          "Not in Universe"               
	1           "Vacation home"                 
	2           "Other residential property"    
	3           "Farm property"                 
	4           "Commercial property"           
	5           "Equipment"                     
	6           "Other"                         
;
label values aritype6 aritypep;
label define aritypep
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values eriat    eriat;  
label define eriat   
	-1          "Not in Universe"               
	1           "Yes"                           
	2           "No"                            
;
label values ariat    ariat;  
label define ariat   
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values eriata   eriata; 
label define eriata  
	-1          "Not in Universe"               
	1           "Yes"                           
	2           "No"                            
;
label values ariata   ariata; 
label define ariata  
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values trimv    trimv;  
label define trimv   
	0           "None or not in universe"       
;
label values arimv    arimv;  
label define arimv   
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values erideb   erideb; 
label define erideb  
	-1          "Not in Universe"               
	1           "Yes"                           
	2           "No"                            
;
label values arideb   arideb; 
label define arideb  
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values tripri   tripri; 
label define tripri  
	0           "None or not in universe"       
;
label values aripri   aripri; 
label define aripri  
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values ertown   ertown; 
label define ertown  
	-1          "Not in Universe"               
	1           "Yes"                           
	2           "No"                            
;
label values artown   artown; 
label define artown  
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values ertnum   ertnum; 
label define ertnum  
	0           "None or not in universe"       
;
label values artnum   artnum; 
label define artnum  
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values erttype1 erttypem;
label define erttypem
	-1          "Not in Universe"               
	1           "Vacation home"                 
	2           "Other residential property"    
	3           "Farm property"                 
	4           "Commercial property"           
	5           "Equipment"                     
	6           "Other"                         
;
label values arttype1 arttypem;
label define arttypem
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values erttype2 erttypek;
label define erttypek
	-1          "Not in Universe"               
	1           "Vacation home"                 
	2           "Other residential property"    
	3           "Farm property"                 
	4           "Commercial property"           
	5           "Equipment"                     
	6           "Other"                         
;
label values arttype2 arttypek;
label define arttypek
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values erttype3 erttypel;
label define erttypel
	-1          "Not in Universe"               
	1           "Vacation home"                 
	2           "Other residential property"    
	3           "Farm property"                 
	4           "Commercial property"           
	5           "Equipment"                     
	6           "Other"                         
;
label values arttype3 arttypel;
label define arttypel
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values erttype4 erttypen;
label define erttypen
	-1          "Not in Universe"               
	1           "Vacation home"                 
	2           "Other residential property"    
	3           "Farm property"                 
	4           "Commercial property"           
	5           "Equipment"                     
	6           "Other"                         
;
label values arttype4 arttypen;
label define arttypen
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values erttype5 erttypeo;
label define erttypeo
	-1          "Not in Universe"               
	1           "Vacation home"                 
	2           "Other residential property"    
	3           "Farm property"                 
	4           "Commercial property"           
	5           "Equipment"                     
	6           "Other"                         
;
label values arttype5 arttypeo;
label define arttypeo
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values erttype6 erttypep;
label define erttypep
	-1          "Not in Universe"               
	1           "Vacation home"                 
	2           "Other residential property"    
	3           "Farm property"                 
	4           "Commercial property"           
	5           "Equipment"                     
	6           "Other"                         
;
label values arttype6 arttypep;
label define arttypep
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values trtmv    trtmv;  
label define trtmv   
	0           "None or not in universe"       
;
label values artmv    artmv;  
label define artmv   
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values ertdeb   ertdeb; 
label define ertdeb  
	-1          "Not in Universe"               
	1           "Yes"                           
	2           "No"                            
;
label values artdeb   artdeb; 
label define artdeb  
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values trtpri   trtpri; 
label define trtpri  
	0           "None or not in universe"       
;
label values artpri   artpri; 
label define artpri  
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values trtsha   trtsha; 
label define trtsha  
	0           "None or not in universe"       
;
label values artsha   artsha; 
label define artsha  
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values tmjp     tmjp;   
label define tmjp    
	0           "None or not in universe"       
;
label values amjp     amjp;   
label define amjp    
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values tmip     tmip;   
label define tmip    
	0           "None or not in universe"       
;
label values amip     amip;   
label define amip    
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values evbunv1  evbunv1l;
label define evbunv1l
	-1          "Not in Universe"               
	1           "In universe"                   
;
label values evbno1   evbno1l;
label define evbno1l 
	-1          "Not in Universe"               
;
label values evbow1   evbow1l;
label define evbow1l 
	0           "Not In Universe"               
;
label values avbow1   avbow1l;
label define avbow1l 
	0           "Not imputed"                   
	1           "Statistical imputed (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values tvbva1   tvbva1l;
label define tvbva1l 
	0           "None or not in universe"       
;
label values avbva1   avbva1l;
label define avbva1l 
	0           "Not imputed"                   
	1           "Statistical imputed (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values tvbde1   tvbde1l;
label define tvbde1l 
	0           "None or not in universe"       
;
label values avbde1   avbde1l;
label define avbde1l 
	0           "Not imputed"                   
	1           "Statistical imputed (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values evbunv2  evbunv2l;
label define evbunv2l
	-1          "Not in Universe"               
	1           "In universe"                   
;
label values evbno2   evbno2l;
label define evbno2l 
	-1          "Not in Universe"               
;
label values evbow2   evbow2l;
label define evbow2l 
	0           "Not In Universe"               
;
label values avbow2   avbow2l;
label define avbow2l 
	0           "Not imputed"                   
	1           "Statistical imputed (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values tvbva2   tvbva2l;
label define tvbva2l 
	0           "None or not in universe"       
;
label values avbva2   avbva2l;
label define avbva2l 
	0           "Not imputed"                   
	1           "Statistical imputed (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values tvbde2   tvbde2l;
label define tvbde2l 
	0           "None or not in universe"       
;
label values avbde2   avbde2l;
label define avbde2l 
	0           "Not imputed"                   
	1           "Statistical imputed (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values emdunv   emdunv; 
label define emdunv  
	-1          "Not in Universe"               
	1           "In universe"                   
;
label values tdonorid tdonorid;
label define tdonorid
	0           "Not in universe or did not"    
	1           "Received data from a donor"    
;
label values ehouspay ehouspay;
label define ehouspay
	-1          "Not in Universe"               
	1           "Yes"                           
	2           "No"                            
;
label values ahouspay ahouspay;
label define ahouspay
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values efoodpay efoodpay;
label define efoodpay
	-1          "Not in Universe"               
	1           "Yes"                           
	2           "No"                            
;
label values afoodpay afoodpay;
label define afoodpay
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values eexppay  eexppay;
label define eexppay 
	-1          "Not in Universe"               
	1           "Yes"                           
	2           "No"                            
;
label values aexppay  aexppay;
label define aexppay 
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values ehhpay   ehhpay; 
label define ehhpay  
	-1          "Not in Universe"               
	1           "Yes"                           
	2           "No"                            
;
label values ahhpay   ahhpay; 
label define ahhpay  
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values ewhopy01 ewhopy0y;
label define ewhopy0y
	-1          "Not in Universe"               
;
label values ewhopy02 ewhopy0k;
label define ewhopy0k
	-1          "Not in Universe"               
;
label values ewhopy03 ewhopy0l;
label define ewhopy0l
	-1          "Not in Universe"               
;
label values ewhopy04 ewhopy0m;
label define ewhopy0m
	-1          "Not in Universe"               
;
label values ewhopy05 ewhopy0n;
label define ewhopy0n
	-1          "Not in Universe"               
;
label values ewhopy06 ewhopy0o;
label define ewhopy0o
	-1          "Not in Universe"               
;
label values ewhopy07 ewhopy0p;
label define ewhopy0p
	-1          "Not in Universe"               
;
label values ewhopy08 ewhopy0q;
label define ewhopy0q
	-1          "Not in Universe"               
;
label values ewhopy09 ewhopy0r;
label define ewhopy0r
	-1          "Not in Universe"               
;
label values ewhopy10 ewhopy1y;
label define ewhopy1y
	-1          "Not in Universe"               
;
label values ewhopy11 ewhopy1k;
label define ewhopy1k
	-1          "Not in Universe"               
;
label values ewhopy12 ewhopy1l;
label define ewhopy1l
	-1          "Not in Universe"               
;
label values ewhopy13 ewhopy1m;
label define ewhopy1m
	-1          "Not in Universe"               
;
label values ewhopy14 ewhopy1n;
label define ewhopy1n
	-1          "Not in Universe"               
;
label values ewhopy15 ewhopy1o;
label define ewhopy1o
	-1          "Not in Universe"               
;
label values ewhopy16 ewhopy1p;
label define ewhopy1p
	-1          "Not in Universe"               
;
label values ewhopy17 ewhopy1q;
label define ewhopy1q
	-1          "Not in Universe"               
;
label values ewhopy18 ewhopy1r;
label define ewhopy1r
	-1          "Not in Universe"               
;
label values ewhopy19 ewhopy1s;
label define ewhopy1s
	-1          "Not in Universe"               
;
label values ewhopy20 ewhopy2y;
label define ewhopy2y
	-1          "Not in Universe"               
;
label values ewhopy21 ewhopy2k;
label define ewhopy2k
	-1          "Not in Universe"               
;
label values ewhopy22 ewhopy2l;
label define ewhopy2l
	-1          "Not in Universe"               
;
label values ewhopy23 ewhopy2m;
label define ewhopy2m
	-1          "Not in Universe"               
;
label values ewhopy24 ewhopy2n;
label define ewhopy2n
	-1          "Not in Universe"               
;
label values ewhopy25 ewhopy2o;
label define ewhopy2o
	-1          "Not in Universe"               
;
label values ewhopy26 ewhopy2p;
label define ewhopy2p
	-1          "Not in Universe"               
;
label values ewhopy27 ewhopy2q;
label define ewhopy2q
	-1          "Not in Universe"               
;
label values ewhopy28 ewhopy2r;
label define ewhopy2r
	-1          "Not in Universe"               
;
label values ewhopy29 ewhopy2s;
label define ewhopy2s
	-1          "Not in Universe"               
;
label values ewhopy30 ewhopy3y;
label define ewhopy3y
	-1          "Not in Universe"               
;
label values awhopy   awhopy; 
label define awhopy  
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values ehltstat ehltstat;
label define ehltstat
	-1          "Not in Universe"               
	1           "Excellent"                     
	2           "Very Good"                     
	3           "Good"                          
	4           "Fair"                          
	5           "Poor"                          
;
label values ahltstat ahltstat;
label define ahltstat
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values ehospsta ehospsta;
label define ehospsta
	-1          "Not in Universe"               
	1           "Yes"                           
	2           "No"                            
;
label values ahospsta ahospsta;
label define ahospsta
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values ehospnit ehospnit;
label define ehospnit
	0           "None or not in universe"       
;
label values ahospnit ahospnit;
label define ahospnit
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values ehreas1  ehreas1l;
label define ehreas1l
	-1          "Not in Universe"               
	1           "Yes"                           
	2           "No"                            
;
label values ahreas1  ahreas1l;
label define ahreas1l
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values ehreas2  ehreas2l;
label define ehreas2l
	-1          "Not in Universe"               
	1           "Yes"                           
	2           "No"                            
;
label values ahreas2  ahreas2l;
label define ahreas2l
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values ehreas3  ehreas3l;
label define ehreas3l
	-1          "Not in Universe"               
	1           "Yes"                           
	2           "No"                            
;
label values ahreas3  ahreas3l;
label define ahreas3l
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values ehreas4  ehreas4l;
label define ehreas4l
	-1          "Not in Universe"               
	1           "Yes"                           
	2           "No"                            
;
label values ahreas4  ahreas4l;
label define ahreas4l
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values ehreas5  ehreas5l;
label define ehreas5l
	-1          "Not in Universe"               
	1           "Yes"                           
	2           "No"                            
;
label values ahreas5  ahreas5l;
label define ahreas5l
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values ehreas6  ehreas6l;
label define ehreas6l
	-1          "Not in Universe"               
	1           "Yes"                           
	2           "No"                            
;
label values ahreas6  ahreas6l;
label define ahreas6l
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values edocnum  edocnum;
label define edocnum 
	0           "None or not in universe"       
;
label values adocnum  adocnum;
label define adocnum 
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values thipay   thipay; 
label define thipay  
	0           "Not in universe or none"       
;
label values ahipay   ahipay; 
label define ahipay  
	0           "Not imputed"                   
	1           "Hot deck"                      
	2           "Hot deck (using unfolding"     
	3           "Logical imputation"            
	4           "Logical imputation (using"     
;
label values epresdrg epresdrg;
label define epresdrg
	-1          "Not in Universe"               
	1           "Yes"                           
	2           "No"                            
;
label values apresdrg apresdrg;
label define apresdrg
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values edalydrg edalydrg;
label define edalydrg
	-1          "Not in Universe"               
	1           "Yes"                           
	2           "No"                            
;
label values adalydrg adalydrg;
label define adalydrg
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values evisdent evisdent;
label define evisdent
	0           "None or not in universe"       
;
label values avisdent avisdent;
label define avisdent
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values edenseal edenseal;
label define edenseal
	-1          "Not in Universe"               
	1           "Yes"                           
	2           "No"                            
;
label values adenseal adenseal;
label define adenseal
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values edis1    edis1l; 
label define edis1l  
	-1          "Not in Universe"               
	1           "Yes"                           
	2           "No"                            
;
label values edis2    edis2l; 
label define edis2l  
	-1          "Not in Universe"               
	1           "Yes"                           
	2           "No"                            
;
label values edis3    edis3l; 
label define edis3l  
	-1          "Not in Universe"               
	1           "Yes"                           
	2           "No"                            
;
label values edis4    edis4l; 
label define edis4l  
	-1          "Not in Universe"               
	1           "Yes"                           
	2           "No"                            
;
label values edis5    edis5l; 
label define edis5l  
	-1          "Not in Universe"               
	1           "Yes"                           
	2           "No"                            
;
label values edis6    edis6l; 
label define edis6l  
	-1          "Not in Universe"               
	1           "Yes"                           
	2           "No"                            
;
label values adis1    adis1l; 
label define adis1l  
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values adis2    adis2l; 
label define adis2l  
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values adis3    adis3l; 
label define adis3l  
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values adis4    adis4l; 
label define adis4l  
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values adis5    adis5l; 
label define adis5l  
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values adis6    adis6l; 
label define adis6l  
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values elostth  elostth;
label define elostth 
	-1          "Not in Universe"               
	1           "Yes"                           
	2           "No"                            
;
label values alostth  alostth;
label define alostth 
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values eallth   eallth; 
label define eallth  
	-1          "Not in Universe"               
	1           "Yes"                           
	2           "No"                            
;
label values aallth   aallth; 
label define aallth  
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values evisdoc  evisdoc;
label define evisdoc 
	0           "None or not in universe"       
;
label values avisdoc  avisdoc;
label define avisdoc 
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values emdspnd  emdspnd;
label define emdspnd 
	-1          "Not in Universe"               
	1           "Yes"                           
	2           "No"                            
;
label values amdspnd  amdspnd;
label define amdspnd 
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values emdspnds emdspnds;
label define emdspnds
	-1          "Not in Universe"               
	1           "Yes"                           
	2           "No"                            
;
label values amdspnds amdspnds;
label define amdspnds
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values edaysick edaysick;
label define edaysick
	0           "None or not in universe"       
;
label values adaysick adaysick;
label define adaysick
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values tmdpay   tmdpay; 
label define tmdpay  
	0           "Not in universe or none"       
;
label values amdpay   amdpay; 
label define amdpay  
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values ereimb   ereimb; 
label define ereimb  
	-1          "Not in Universe"               
	1           "Total actual Cost"             
	2           "Got Reimbursed"                
	3           "Expects to get reimbursed but has"
;
label values areimb   areimb; 
label define areimb  
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values treimbur treimbur;
label define treimbur
	0           "None or not in universe"       
;
label values areimbur areimbur;
label define areimbur
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values ehspstas ehspstas;
label define ehspstas
	-1          "Not in Universe"               
	1           "Yes"                           
	2           "No"                            
;
label values ahspstas ahspstas;
label define ahspstas
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values eprsdrgs eprsdrgs;
label define eprsdrgs
	-1          "Not in Universe"               
	1           "Yes"                           
	2           "No"                            
;
label values aprsdrgs aprsdrgs;
label define aprsdrgs
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values evsdents evsdents;
label define evsdents
	-1          "Not in Universe"               
	1           "Yes"                           
	2           "No"                            
;
label values avsdents avsdents;
label define avsdents
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values evsdocs  evsdocs;
label define evsdocs 
	-1          "Not in Universe"               
	1           "Yes"                           
	2           "No"                            
;
label values avsdocs  avsdocs;
label define avsdocs 
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values enowkyr  enowkyr;
label define enowkyr 
	-1          "Not in Universe"               
	1           "A year or longer"              
	2           "less than a year"              
;
label values anowkyr  anowkyr;
label define anowkyr 
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values ewkfutr  ewkfutr;
label define ewkfutr 
	-1          "Not in Universe"               
	1           "Yes"                           
	2           "No"                            
;
label values awkfutr  awkfutr;
label define awkfutr 
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values trmoops  trmoops;
label define trmoops 
	0           "None or not in universe"       
;
label values enoindnt enoindnt;
label define enoindnt
	-1          "Not in Universe"               
	1           "Yes"                           
	2           "No"                            
;
label values anoindnt anoindnt;
label define anoindnt
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values enoindoc enoindoc;
label define enoindoc
	-1          "Not in Universe"               
	1           "Yes"                           
	2           "No"                            
;
label values anoindoc anoindoc;
label define anoindoc
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values enointrt enointrt;
label define enointrt
	-1          "Not in Universe"               
	1           "Yes"                           
	2           "No"                            
;
label values anointrt anointrt;
label define anointrt
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values enoinchk enoinchk;
label define enoinchk
	-1          "Not in Universe"               
	1           "Yes"                           
	2           "No"                            
;
label values anoinchk anoinchk;
label define anoinchk
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values enoindrg enoindrg;
label define enoindrg
	-1          "Not in Universe"               
	1           "Yes"                           
	2           "No"                            
;
label values anoindrg anoindrg;
label define anoindrg
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values enoinpay enoinpay;
label define enoinpay
	-1          "Not in Universe"               
	1           "Free"                          
	2           "Paid something"                
	3           "Both (if respondent volunteers)"
;
label values anoinpay anoinpay;
label define anoinpay
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values enoindis enoindis;
label define enoindis
	-1          "Not in Universe"               
	1           "Full price"                    
	2           "Reduced price"                 
	3           "Don't know"                    
;
label values anoindis anoindis;
label define anoindis
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values enoininc enoininc;
label define enoininc
	-1          "Not in Universe"               
	1           "Yes"                           
	2           "No"                            
;
label values anoininc anoininc;
label define anoininc
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values enoincln enoincln;
label define enoincln
	-1          "Not in Universe"               
	1           "Yes"                           
	2           "No"                            
;
label values enoiner  enoiner;
label define enoiner 
	-1          "Not in Universe"               
	1           "Yes"                           
	2           "No"                            
;
label values enoinhsp enoinhsp;
label define enoinhsp
	-1          "Not in Universe"               
	1           "Yes"                           
	2           "No"                            
;
label values enoinva  enoinva;
label define enoinva 
	-1          "Not in Universe"               
	1           "Yes"                           
	2           "No"                            
;
label values enoindr  enoindr;
label define enoindr 
	-1          "Not in Universe"               
	1           "Yes"                           
	2           "No"                            
;
label values enoindds enoindds;
label define enoindds
	-1          "Not in Universe"               
	1           "Yes"                           
	2           "No"                            
;
label values enoinoth enoinoth;
label define enoinoth
	-1          "Not in Universe"               
	1           "Yes"                           
	2           "No"                            
;
label values anoinloc anoinloc;
label define anoinloc
	0           "Not imputed"                   
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values eapvunv  eapvunv;
label define eapvunv 
	-1          "Not in Universe"               
	1           "In universe"                   
;
label values epvwk1   epvwk1l;
label define epvwk1l 
	-1          "Not in Universe"               
	1           "Yes"                           
	2           "No"                            
;
label values epvwk2   epvwk2l;
label define epvwk2l 
	-1          "Not in Universe"               
	1           "Yes"                           
	2           "No"                            
;
label values epvwk3   epvwk3l;
label define epvwk3l 
	-1          "Not in Universe"               
	1           "Yes"                           
	2           "No"                            
;
label values epvwk4   epvwk4l;
label define epvwk4l 
	-1          "Not in Universe"               
	1           "Yes"                           
	2           "No"                            
;
label values epvwk5   epvwk5l;
label define epvwk5l 
	-1          "Not in Universe"               
	1           "Yes"                           
	2           "No"                            
;
label values apvwk    apvwk;  
label define apvwk   
	0           "No imputation"                 
	1           "Statistical imputation (hot deck)"
	2           "Cold deck"                     
	3           "Logical imputation (derivation)"
;
label values epvmilwk epvmilwk;
label define epvmilwk
	-1          "Not in Universe"               
;
label values apvmilwk apvmilwk;
label define apvmilwk
	0           "No imputation"                 
	1           "Statistical imputation (hot deck)"
	2           "Cold deck"                     
	3           "Logical imputation (derivation)"
;
label values epvpaprk epvpaprk;
label define epvpaprk
	-1          "Not in Universe"               
	1           "Yes"                           
	2           "No"                            
;
label values apvpaprk apvpaprk;
label define apvpaprk
	0           "No imputation"                 
	1           "Statistical imputation (hot deck)"
	2           "Cold deck"                     
	3           "Logical imputation (derivation)"
;
label values epvpaywk epvpaywk;
label define epvpaywk
	0           "Not In Universe"               
;
label values apvpaywk apvpaywk;
label define apvpaywk
	0           "No imputation"                 
	1           "Statistical imputation (hot deck)"
	2           "Cold deck"                     
	3           "Logical imputation (derivation)"
;
label values epvcomut epvcomut;
label define epvcomut
	0           "Not In Universe"               
;
label values apvcomut apvcomut;
label define apvcomut
	0           "No imputation"                 
	1           "Statistical imputation (hot deck)"
	2           "Cold deck"                     
	3           "Logical imputation (derivation)"
;
label values epvwkexp epvwkexp;
label define epvwkexp
	-1          "Not in Universe"               
	1           "Yes"                           
	2           "No"                            
;
label values apvwkexp apvwkexp;
label define apvwkexp
	0           "No imputation"                 
	1           "Statistical imputation (hot deck)"
	2           "Cold deck"                     
	3           "Logical imputation (derivation)"
;
label values epvanexp epvanexp;
label define epvanexp
	0           "Not In Universe"               
;
label values apvanexp apvanexp;
label define apvanexp
	0           "No imputation"                 
	1           "Statistical imputation (hot deck)"
	2           "Cold deck"                     
	3           "Logical imputation (derivation)"
;
label values epvchild epvchild;
label define epvchild
	-1          "Not in Universe"               
	1           "Yes"                           
	2           "No"                            
;
label values apvchild apvchild;
label define apvchild
	0           "No imputation"                 
	1           "Statistical imputation (hot deck)"
	2           "Cold deck"                     
	3           "Logical imputation (derivation)"
;
label values epvmancd epvmancd;
label define epvmancd
	-1          "Not in Universe"               
;
label values apvmancd apvmancd;
label define apvmancd
	0           "No imputation"                 
	1           "Statistical imputation (hot deck)"
	2           "Cold deck"                     
	3           "Logical imputation (derivation)"
;
label values epvmosup epvmosup;
label define epvmosup
	-1          "Not in Universe"               
	1           "Yes"                           
	2           "No"                            
;
label values apvmosup apvmosup;
label define apvmosup
	0           "No imputation"                 
	1           "Statistical imputation (hot deck)"
	2           "Cold deck"                     
	3           "Logical imputation (derivation)"
;
label values tpvchpa1 tpvchpap;
label define tpvchpap
	0           "None or not in universe"       
;
label values tpvchpa2 tpvchpak;
label define tpvchpak
	0           "None or not in universe"       
;
label values tpvchpa3 tpvchpal;
label define tpvchpal
	0           "None or not in universe"       
;
label values tpvchpa4 tpvchpam;
label define tpvchpam
	0           "None or not in universe"       
;
label values apvchpa  apvchpa;
label define apvchpa 
	0           "No imputation"                 
	1           "Statistical imputation (hot deck)"
	2           "Cold deck"                     
	3           "Logical imputation (derivation)"
;
label values epvccarr epvccarr;
label define epvccarr
	-1          "Not in Universe"               
	1           "Yes"                           
	2           "No"                            
;
label values apvccarr apvccarr;
label define apvccarr
	0           "No imputation"                 
	1           "Statistical imputation (hot deck)"
	2           "Cold deck"                     
	3           "Logical imputation (derivation)"
;
label values tpvccfp1 tpvccfpr;
label define tpvccfpr
	0           "None or not in universe"       
;
label values apvccfp1 apvccfpr;
label define apvccfpr
	0           "No imputation"                 
	1           "Statistical imputation (hot deck)"
	2           "Cold deck"                     
	3           "Logical imputation (derivation)"
;
label values tpvccfp2 tpvccfpk;
label define tpvccfpk
	0           "None or not in universe"       
;
label values apvccfp2 apvccfpk;
label define apvccfpk
	0           "No imputation"                 
	1           "Statistical imputation (hot deck)"
	2           "Cold deck"                     
	3           "Logical imputation (derivation)"
;
label values tpvccfp3 tpvccfpl;
label define tpvccfpl
	0           "None or not in universe"       
;
label values apvccfp3 apvccfpl;
label define apvccfpl
	0           "No imputation"                 
	1           "Statistical imputation (hot deck)"
	2           "Cold deck"                     
	3           "Logical imputation (derivation)"
;
label values tpvccfp4 tpvccfpm;
label define tpvccfpm
	0           "None or not in universe"       
;
label values apvccfp4 apvccfpm;
label define apvccfpm
	0           "No imputation"                 
	1           "Statistical imputation (hot deck)"
	2           "Cold deck"                     
	3           "Logical imputation (derivation)"
;
label values epvccoth epvccoth;
label define epvccoth
	-1          "Not in Universe"               
	1           "Yes"                           
	2           "No"                            
;
label values apvccoth apvccoth;
label define apvccoth
	0           "No imputation"                 
	1           "Statistical imputation (hot deck)"
	2           "Cold deck"                     
	3           "Logical imputation (derivation)"
;
label values epvcwho1 epvcwhoh;
label define epvcwhoh
	-1          "Not in Universe"               
	1           "Yes"                           
	2           "No"                            
;
label values epvcwho2 epvcwhok;
label define epvcwhok
	-1          "Not in Universe"               
	1           "Yes"                           
	2           "No"                            
;
label values epvcwho3 epvcwhol;
label define epvcwhol
	-1          "Not in Universe"               
	1           "Yes"                           
	2           "No"                            
;
label values epvcwho4 epvcwhom;
label define epvcwhom
	-1          "Not in Universe"               
	1           "Yes"                           
	2           "No"                            
;
label values epvcwho5 epvcwhon;
label define epvcwhon
	-1          "Not in Universe"               
	1           "Yes"                           
	2           "No"                            
;
label values apvcwho  apvcwho;
label define apvcwho 
	0           "No imputation"                 
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;
label values epvdays  epvdays;
label define epvdays 
	-1          "Not in Universe"               
;
label values epvweeks epvweeks;
label define epvweeks
	-1          "Not in Universe"               
;
label values epvmnths epvmnths;
label define epvmnths
	-1          "Not in Universe"               
;
label values apvdwm   apvdwm; 
label define apvdwm  
	0           "No imputation"                 
	1           "Statistical imputation (hot deck)"
	2           "Cold deck imputation"          
	3           "Logical imputation (derivation)"
;

#delimit cr
desc,short

****************************************************************
*to make consistent with CEPR files:
sort ssuid shhadid eentaid epppnum swave srotaton 

rename spanel panel
notes panel: SIPP: spanel (96,01,04)
* gen unique id to merge
sort ssuid eentaid epppnum
egen id = concat(ssuid eentaid epppnum)


rename swave wave
rename srotaton rot
sort id rot wave

keep id rot wave ssuid eentaid epppnum   thomeamt thh* tal*  tmhpr
replace thhdebt  = thhdebt  + tmhpr
replace thhscdbt = thhscdbt + tmhpr
// note: should be the case thhdebt = thhscdbt + thhuscbt? 

saveold "${save_dir}sippp08putm7.dta", replace


** Copyright 2012 shared by the National Bureau of Economic Research and Jean Roth ;

** National Bureau of Economic Research. ;
** 1050 Massachusetts Avenue ;
** Cambridge, MA 02138 ;
** jroth@nber.org ;

** This program and all programs referenced in it are free software. You ;
** can redistribute the program or modify it under the terms of the GNU ;
** General Public License as published by the Free Software Foundation; 
** either version 2 of the License, or (at your option) any later version. ;

** This program is distributed in the hope that it will be useful, ;
** but WITHOUT ANY WARRANTY, without even the implied warranty of ;
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the ;
** GNU General Public License for more details. ;

** You should have received a copy of the GNU General Public License ;
** along with this program, if not, write to the Free Software ;
** Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA. ;
