

universe        = vanilla
getenv          = true
executable      = runLimit.sh
output          = job.Limit.$(Process).out
error           = job.Limit.$(Process).err
log             = job.Limit.$(Process).log
priority        = 10

Notification    = Error

Requirements =  ( SlotID == 1 || SlotID == 2 || SlotId == 3 || SlotId == 4 || SlotId5== 5) 



should_transfer_files = YES
when_to_transfer_output = ON_EXIT
transfer_input_files = Histfitter.tgz
transfer_output_files = 

Arguments = $(process)

queue 1