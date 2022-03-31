#include "AtParsers.h"

ClassImp(AtTools::AtParsers)

AtTools::AtParsers::AtParsers()
{

}

AtTools::AtParsers::~AtParsers()
{

}

Int_t AtTools::AtParsers::ParseIonFitXML(TString filename)
{

  std::shared_ptr<TDOMParser> parser =std::make_shared<TDOMParser>();
  
  Int_t parsecode = parser->ParseFile(filename);
 
      if (parsecode < 0) {
	std::cerr << parser->GetParseCodeMessage(parsecode) << "\n";
         return -1;
      }

      TXMLNode * node = parser->GetXMLDocument()->GetRootNode();
      ParseIonList(node);
      return 0;    
} 

void AtTools::AtParsers::ParseIonList(TXMLNode *node)
{
 
       for (; node; node = node->GetNextNode()) {

	 if (node->GetNodeType() == TXMLNode::kXMLElementNode) { // Element Node
	   
	   if (strcmp(node->GetNodeName(), "ion") == 0) {
	   
	     Int_t id=0;
               if (node->HasAttributes()) {
                  TList *attrList = node->GetAttributes();
                  TXMLAttr *attr = 0;
                  TIter next(attrList);
                  while ((attr=(TXMLAttr*)next())) {
                     if (strcmp(attr->GetName(), "ID") == 0) {
                        id = atoi(attr->GetValue());			
                        break;
                     }
                  }
               }//Attributes

	       ionList.push_back(ParseIon(node->GetChildren(), id));
	       
	   }//ion
	   
         }//Element node
	 
	 ParseIonList(node->GetChildren());
      }

}  

AtTools::IonFitInfo AtTools::AtParsers::ParseIon(TXMLNode *node, Int_t id)
{

   std::string ionName;
   UInt_t PDG;
   Double_t mass;
   Int_t atomicNumber;
   UInt_t MassNumber;
   std::string eLossFile;

  
      for ( ; node; node = node->GetNextNode()) {
         if (node->GetNodeType() == TXMLNode::kXMLElementNode) { // Element Node
	   if (strcmp(node->GetNodeName(), "Ion") == 0){
               ionName = node->GetText();
	       
	   }
            if (strcmp(node->GetNodeName(), "PDG") == 0){
	      PDG = std::atoi(node->GetText());
	      
	    }
	    if (strcmp(node->GetNodeName(), "Mass") == 0){
	      mass = std::stod(node->GetText());
	      
            }
            if (strcmp(node->GetNodeName(), "AtomicNumber") == 0){
	      atomicNumber = std::atoi(node->GetText());
	      
            }
            if (strcmp(node->GetNodeName(), "MassNumber") == 0){
	      MassNumber = std::atoi(node->GetText());
	      
            }
            if (strcmp(node->GetNodeName(), "ElossFile") == 0){
	      eLossFile = node->GetText();
	      
            }

	 }
     }
 
      AtTools::IonFitInfo ionInfo={ionName, PDG, mass,atomicNumber,MassNumber,eLossFile};
         
      return ionInfo;
}  
