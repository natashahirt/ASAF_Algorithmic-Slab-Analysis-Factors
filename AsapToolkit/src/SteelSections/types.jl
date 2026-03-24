abstract type AbstractSection end
abstract type TorsionAllowed <: AbstractSection end

"""
Wide-flange sections
"""
@SteelSections W TorsionAllowed Wfields Wrange

"""
Channels
"""
@SteelSections C TorsionAllowed Cfields Crange

"""
Angles
"""
@SteelSections L TorsionAllowed Lfields Lrange

"""
Double angles
"""
@SteelSections LL AbstractSection LLfields LLrange

"""
WT sections
"""
@SteelSections WT TorsionAllowed Wfields WTrange

"""
HSS rectangular
"""
@SteelSections HSSRect TorsionAllowed HSSRectfields HSSRectrange

"""
HSS round
"""
@SteelSections HSSRound TorsionAllowed HSSRoundfields HSSRoundrange
