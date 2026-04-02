#include "AtFitter.h"

#include "AtPatternEvent.h"
#include "AtTrackingEvent.h"

void EventFit::AtFitter::FitEvent(AtTrackingEvent *trackingEvent, AtPatternEvent *patternEvent,
                                  AtFitMetadata *fitMetadata, AtRawEvent *rawEvent, AtEvent *event)
{
   // Check for nullptr.
   if (trackingEvent == nullptr) {
      LOG(error) << " Tracking event is nullptr! The fitter can not fit this event. Maybe the tracking event is not "
                    "being constructed properly in the fitter task.";
      return;
   }

   if (patternEvent == nullptr) {
      LOG(error) << " Pattern event is nullptr! The fitter can not fit this event.";
      return;
   }

   // Extract the candidate AtTracks. If there are not any tracks, return earlier.
   LOG(info) << "FitEvent: copying tracks...";
   std::vector<AtTrack> tracks = patternEvent->GetTrackCand();
   if (!tracks.size())
      return;
   LOG(info) << "FitEvent: " << tracks.size() << " tracks copied, saving to tracking event...";

   // Save the original AtTracks to the AtTrackingEvent.
   trackingEvent->SetTrackArray(&tracks);
   LOG(info) << "FitEvent: fitting tracks...";

   // Iterate over the AtTracks and store the AtFittedTracks in the AtTrackingEvent.
   for (size_t i = 0; i < tracks.size(); ++i) {
      LOG(info) << "FitEvent: fitting track " << i << " with " << tracks[i].GetHitClusterArray()->size() << " clusters";
      std::unique_ptr<AtFittedTrack> fittedTrack(GetFittedTrack(&tracks[i], fitMetadata, rawEvent, event));
      if (fittedTrack) {
         LOG(info) << "FitEvent: track " << i << " fitted successfully";
         trackingEvent->AddFittedTrack(std::move(fittedTrack));
      } else {
         LOG(info) << "FitEvent: track " << i << " returned null";
      }
   }
}
